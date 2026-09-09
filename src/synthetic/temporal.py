"""Deterministic temporal forecasting benchmarks; no biological claims."""
import numpy as np
def stationary(n_patients=30,seed=0,heavy_tail=False):
    r=np.random.default_rng(seed); A=np.array([[.7,.1],[0,.6]]); b=np.array([.2,-.1]); z=r.normal(size=(n_patients,2)); e=r.normal(size=(n_patients,2))*.1
    if heavy_tail: e=r.standard_t(3,size=(n_patients,2))*.1
    return z,z@A.T+b+e,np.array(['A']*n_patients),np.ones(n_patients),np.array([f'p{i}' for i in range(n_patients)])
def irregular(seed=1):
    return {'patient_id':['p','p','p'],'time':np.array([0.,1.,3.]),'state':[np.array([0.]),np.array([1.]),np.array([3.])],'treatment':['A','A','A']}
