import numpy as np
import pandas as pd
import pytest
from src.validation.temporal import make_forecast_example,HistoricalScaler
from src.representation.encoder import VersionedEncoder
from src.dynamics.discrete import ProbabilisticTransition
from src.synthetic.temporal import stationary

def frame():
 return pd.DataFrame({'patient_id':['A']*4+['B']*4,'time':[0,1,2,3]*2,'state':[np.array([i]) for i in [0,1,2,3,0,1,2,3]],'treatment':['A']*8})
def test_history_target_contract():
 e=make_forecast_example(frame(),'A',2); assert e.history.time.max()<e.prediction_time; assert e.target.time==2; assert e.delta_t==1
 assert 2 not in e.history.time.tolist() and not (e.history.time>2).any()
def test_target_and_later_rows_cannot_enter_history():
 with pytest.raises(ValueError): make_forecast_example(frame(),'A',0)
def test_future_normalization_does_not_change_historical_parameters():
 historical=np.array([[0.],[1.],[2.]]); future=np.array([[1000.]])
 a=HistoricalScaler().fit(historical); b=HistoricalScaler().fit(np.vstack([historical,future])); assert np.array_equal(a.mean_,[1.]); assert not np.array_equal(a.mean_,b.mean_)
 c=HistoricalScaler().fit(historical); c.transform(future); assert np.array_equal(a.mean_,c.mean_)
def test_encoder_is_frozen_against_future_rows():
 e=VersionedEncoder.from_yaml('configs/representation/latent_state_v1.yaml'); h=e.encoder_hash
 historical={'MKI67':1.,'TOP2A':2.,'UBE2C':3.}; future={'MKI67':1e12,'TOP2A':-1e12,'UBE2C':1e12}
 before=e.encode(historical,sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=100,available_modalities=['scRNA-seq'])
 _=e.encode(future,sample_id='future',patient_id='p',cohort_id='c',timepoint='t+1',n_cells=100,available_modalities=['scRNA-seq'])
 after=e.encode(historical,sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=100,available_modalities=['scRNA-seq'])
 assert e.encoder_hash==h
 for name in before.features:
  assert before.features[name].n_cells==after.features[name].n_cells
  assert (np.isnan(before.features[name].value) and np.isnan(after.features[name].value)) or before.features[name].value==after.features[name].value
 e2=VersionedEncoder.from_yaml('configs/representation/latent_state_v1.yaml'); assert e2.encoder_hash==h
def test_future_transitions_do_not_change_historical_covariance():
 z,y,t,dt,p=stationary(12); m=ProbabilisticTransition().fit(z,y,t,dt,p); cov=np.asarray(m.residual_covariance['A'])
 z2,y2,t2,dt2,p2=stationary(12,99); p2=np.array([f'f{i}' for i in range(12)]); m2=ProbabilisticTransition().fit(np.vstack([z,z2]),np.vstack([y,y2]),np.r_[t,t2],np.r_[dt,dt2],np.r_[p,p2]); assert np.allclose(cov,np.asarray(m.residual_covariance['A'])) and not np.allclose(cov,np.asarray(m2.residual_covariance['A']))
