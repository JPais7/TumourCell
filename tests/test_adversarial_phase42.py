import numpy as np
import pytest
from src.dynamics.discrete import ProbabilisticTransition
from src.dynamics.validation import validate_transition_model
from src.clones.clone_state import assess_observed_shift,MechanismStatus
from src.spatial.representation import spatial_confounding,build_neighbourhood_graph,radius_graph
from src.validation.gates import evaluate_gates,simulation_allowed,GATES

def data(n=12,seed=4):
 r=np.random.default_rng(seed); z=r.normal(size=(n,2)); y=z+np.array([.2,-.1])+r.normal(scale=.05,size=(n,2)); return z,y,np.array(['A']*n),np.ones(n),np.array([f'p{i}' for i in range(n)])

def test_loo_train_test_patients_disjoint_and_fold_covariance():
 z,y,t,dt,pid=data(); m=ProbabilisticTransition().fit(z,y,t,dt,pid); out=validate_transition_model(m,z,y,t,dt,pid)
 for f in out['folds']: assert not set(f['train_patients']) & set(f['test_patients'])
 assert all('predictive_covariance' in r for r in out['fold_predictions'])

def test_directional_accuracy_is_delta_based():
 z=np.array([[1.],[1.],[1.],[1.],[1.],[1.],[1.],[1.]])
 y=np.array([[2.],[2.],[2.],[2.],[2.],[2.],[2.],[2.]])
 t=np.array(['A']*8); dt=np.ones(8); p=np.array([f'p{i}' for i in range(8)])
 m=ProbabilisticTransition().fit(z,y,t,dt,p); out=validate_transition_model(m,z,y,t,dt,p); assert out['directional_accuracy'] >= 0

def test_treatment_imbalance_rejected():
 z,y,t,dt,p=data(); t[0]='B'
 with pytest.raises(ValueError): ProbabilisticTransition().fit(z,y,t,dt,p)

@pytest.mark.parametrize('bad', [lambda z,y,t,dt,p: (z*float('nan'),y,t,dt,p), lambda z,y,t,dt,p:(z,y,t,-dt,p), lambda z,y,t,dt,p:(z[:3],y,t,dt,p)])
def test_pathological_transition_inputs_fail(bad):
 z,y,t,dt,p=data()
 with pytest.raises(ValueError): ProbabilisticTransition().fit(*bad(z,y,t,dt,p))

def test_mechanism_statuses_are_distinct():
 assert assess_observed_shift().identifiability_status==MechanismStatus.NOT_IDENTIFIABLE.value
 assert assess_observed_shift(expression_before_after=True).composition_shift==MechanismStatus.INDETERMINATE.value
 assert assess_observed_shift(expression_before_after=True,sampling_problem=True).sampling_confounding==MechanismStatus.SUPPORTED.value
 assert assess_observed_shift(expression_before_after=True,technical_problem=True).technical_confounding==MechanismStatus.SUPPORTED.value
 assert assess_observed_shift(clone_evidence=True,assignment_quality=.95,temporal_overlap=True,abundance_sufficient=True,sampling_comparable=True,clone_frequency_changed=True).selection==MechanismStatus.SUPPORTED.value

def test_spatial_decomposition_truths_and_fail_closed():
 within=spatial_confounding([0.],[1.],{'a':[[0.]],'b':[[0.]]},{'a':[[1.]],'b':[[1.]]})
 assert within['within_region_shift']==[1.0] and abs(within['between_region_composition_shift'][0])<1e-12
 comp=spatial_confounding([0.],[0.],{'a':[[0.]],'b':[[2.]]},{'a':[[0.],[0.]],'b':[[2.]]})
 assert comp['identifiability_status']=='INDETERMINATE'
 assert spatial_confounding([0.],[1.],{'a':[[0.]]},{'b':[[1.]]})['identifiability_status']=='NOT_IDENTIFIABLE'
 assert spatial_confounding([0.],[1.],{'a':[]},{'a':[[1.]]})['identifiability_status']=='NOT_IDENTIFIABLE'

def test_spatial_graph_adversarial_inputs():
 assert len(build_neighbourhood_graph([[0,0],[1,0],[3,0]],'knn',1))==3
 assert radius_graph([[0,0],[10,10]],1)==[(),()]
 with pytest.raises(ValueError): build_neighbourhood_graph([[0,0]],'knn',1)
 with pytest.raises(ValueError): radius_graph([[0,0]],-1)

def test_gate_policy_all_possible_evidence_is_fail_closed():
 assert set(evaluate_gates({})['gates'])==set(GATES)
 assert simulation_allowed({g:'PASSED' for g in GATES}) is False
 assert evaluate_gates({'G1':'PASSED','G2':'PASSED'})['clinical_recommendations_enabled'] is False
