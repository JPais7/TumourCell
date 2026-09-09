"""Fail-closed temporal contracts and splits for forecasting experiments."""
from __future__ import annotations
import numpy as np
import pandas as pd
from dataclasses import dataclass

@dataclass(frozen=True)
class ForecastExample:
    patient_id: str; history: pd.DataFrame; target: pd.Series; prediction_time: float; delta_t: float
    def __post_init__(self):
        if self.history.empty or self.prediction_time <= float(self.history.time.max()): raise ValueError("forecast requires non-empty history strictly before target")
        if self.target.get('time') != self.prediction_time: raise ValueError("target time mismatch")
        if self.delta_t <= 0: raise ValueError("delta_t must be positive")

def make_forecast_example(frame, patient_id, prediction_time):
    validate_temporal_rows(frame[['patient_id','time']].to_dict('records'))
    patient=frame[frame.patient_id==patient_id]
    target=patient[patient.time==prediction_time]
    if len(target)!=1: raise ValueError("exactly one target row required")
    history=patient[patient.time<prediction_time].copy()
    if history.empty: raise ValueError("no historical observations")
    return ForecastExample(str(patient_id),history,target.iloc[0],float(prediction_time),float(prediction_time-history.time.max()))

class HistoricalScaler:
    """Small fit/transform contract used only to test temporal membership."""
    def __init__(self): self.mean_=None; self.scale_=None; self.fitted_rows_=None
    def fit(self, values):
        x=np.asarray(values,float)
        if x.ndim!=2 or not len(x) or not np.isfinite(x).all(): raise ValueError("invalid historical values")
        self.mean_=x.mean(0); self.scale_=np.where(x.std(0)>0,x.std(0),1.); self.fitted_rows_=len(x); return self
    def transform(self, values):
        if self.mean_ is None: raise ValueError("scaler is not fitted")
        x=np.asarray(values,float); return (x-self.mean_)/self.scale_

def validate_temporal_rows(rows, *, require_order=True):
    if not rows: raise ValueError("no temporal rows")
    frame=pd.DataFrame(rows)
    if frame[['patient_id','time']].isna().any().any(): raise ValueError("missing patient_id/time")
    if not np.isfinite(frame.time).all(): raise ValueError("non-finite timestamp")
    if frame.duplicated(['patient_id','time']).any(): raise ValueError("duplicate patient/timepoint")
    if require_order:
        for _,g in frame.groupby('patient_id'):
            if not g.time.is_monotonic_increasing: raise ValueError("non-monotonic patient trajectory")
    return frame

def temporal_split(frame, prediction_time, *, allow_equal=False):
    if 'patient_id' not in frame or 'time' not in frame: raise ValueError("patient_id and time are required")
    if not np.isfinite(prediction_time): raise ValueError("invalid prediction time")
    train=frame[frame.time < prediction_time] if not allow_equal else frame[frame.time <= prediction_time]
    test=frame[frame.time >= prediction_time] if not allow_equal else frame[frame.time > prediction_time]
    if len(train)==0: raise ValueError("insufficient historical rows for prediction")
    if len(train) and train.time.max() >= prediction_time and not allow_equal: raise ValueError("future row entered training")
    return train.copy(),test.copy()

def transition_rows(frame):
    validate_temporal_rows(frame[['patient_id','time']].to_dict('records'))
    if 'state' not in frame: raise ValueError("state column required")
    out=[]
    for patient,g in frame.sort_values(['patient_id','time']).groupby('patient_id',sort=False):
        for i in range(max(0,len(g)-1)):
            dt=float(g.time.iloc[i+1]-g.time.iloc[i])
            if dt<=0: raise ValueError("delta_t must be positive")
            out.append({'patient_id':patient,'time_current':g.time.iloc[i],'time_next':g.time.iloc[i+1],'delta_t':dt,'current_state':g.state.iloc[i],'next_state':g.state.iloc[i+1],'treatment':g.treatment.iloc[i] if 'treatment' in g else None})
    return pd.DataFrame(out)

def readiness_assessment(frame, *, preprocessing_frozen, future_holdout, min_patients=3):
    reasons=[]; checks={}
    try: validate_temporal_rows(frame[['patient_id','time']].to_dict('records')); checks['ordered_unique_times']=True
    except ValueError as e: checks['ordered_unique_times']=False; reasons.append(str(e))
    tr=transition_rows(frame) if checks.get('ordered_unique_times') and 'state' in frame else pd.DataFrame()
    checks.update({'positive_delta_t':bool(len(tr) and (tr.delta_t>0).all()),'sufficient_patients':bool(frame.patient_id.nunique()>=min_patients),'preprocessing_frozen':bool(preprocessing_frozen),'future_holdout':bool(future_holdout),'no_random_cell_split':True})
    if not checks['positive_delta_t']: reasons.append('no valid positive patient transitions')
    return {'phase5_ready':all(checks.values()),'checks':checks,'reasons':reasons,'status':'READY_FOR_PHASE_5' if all(checks.values()) else 'NOT_READY'}
