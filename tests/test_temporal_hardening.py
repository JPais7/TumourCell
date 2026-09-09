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
 e2=VersionedEncoder.from_yaml('configs/representation/latent_state_v1.yaml'); assert e2.encoder_hash==h
def test_future_transitions_do_not_change_historical_covariance():
 z,y,t,dt,p=stationary(12); m=ProbabilisticTransition().fit(z,y,t,dt,p); cov=np.asarray(m.residual_covariance['A'])
 z2,y2,t2,dt2,p2=stationary(12,99); p2=np.array([f'f{i}' for i in range(12)]); m2=ProbabilisticTransition().fit(np.vstack([z,z2]),np.vstack([y,y2]),np.r_[t,t2],np.r_[dt,dt2],np.r_[p,p2]); assert np.allclose(cov,np.asarray(m.residual_covariance['A'])) and not np.allclose(cov,np.asarray(m2.residual_covariance['A']))
