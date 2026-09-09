import numpy as np
import pandas as pd
import pytest
from src.validation.temporal import validate_temporal_rows,temporal_split,transition_rows,readiness_assessment
from src.synthetic.temporal import stationary,irregular
from src.dynamics.discrete import ProbabilisticTransition

def test_strict_temporal_split_excludes_future():
 f=pd.DataFrame({'patient_id':['p','p','p'],'time':[0.,1.,2.]}); tr,te=temporal_split(f,2.)
 assert tr.time.max()<2 and te.time.min()>=2
 with pytest.raises(ValueError): temporal_split(f,-1.)
def test_temporal_malformed_rows_fail():
 for rows in [[],[{'patient_id':'p','time':None}],[{'patient_id':'p','time':1},{'patient_id':'p','time':1}],[{'patient_id':'p','time':2},{'patient_id':'p','time':1}]]:
  with pytest.raises(ValueError): validate_temporal_rows(rows)
def test_irregular_delta_is_preserved():
 d=irregular(); out=transition_rows(pd.DataFrame(d)); assert out.delta_t.tolist()==[1.,2.]
def test_synthetic_forecasting_has_known_shape_and_no_real_claim():
 z,y,t,dt,p=stationary(); m=ProbabilisticTransition().fit(z,y,t,dt,p); assert m.predict_mean(z[0],'A',1).shape==(1,2)
def test_future_rows_do_not_enter_earlier_split():
 f=pd.DataFrame({'patient_id':['a','a','b','b'],'time':[0.,2.,0.,2.]}); tr,te=temporal_split(f,2.)
 assert set(tr.patient_id)=={'a','b'} and all(tr.time<2)
def test_readiness_fails_without_future_holdout():
 f=pd.DataFrame({'patient_id':['a','a','b','b','c','c'],'time':[0,1,0,1,0,1],'state':[np.array([0]),np.array([1])]*3})
 assert readiness_assessment(f,preprocessing_frozen=True,future_holdout=False)['phase5_ready'] is False
def test_duplicate_patient_transition_rejected():
 z=np.zeros((6,1)); y=np.ones((6,1)); t=['A']*6; dt=[1]*6
 with pytest.raises(ValueError): ProbabilisticTransition().fit(z,y,t,dt,['p','p','q','r','s','u'])
