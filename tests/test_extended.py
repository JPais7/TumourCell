import hashlib,tempfile,unittest
from pathlib import Path
import numpy as np
from src.clones.clone_assignment import assign_clone
from src.dynamics.inference import gaussian_posterior
from src.experiments import sha256
from src.observation.model import LinearGaussianObservation
from src.perturbations import Perturbation
from src.representation.latent_state import TumourStateDistribution
from src.spatial.representation import radius_graph,neighbourhood_context,validate_coordinates
from src.validation.adjustment import three_way_adjustment
from src.validation.patient_level import cell_threshold_sensitivity
from src.validation.permutation import negative_controls,permutation_test
from src.dynamics.discrete import ProbabilisticTransition
from src.dynamics.validation import validate_transition_model
from src.spatial.representation import neighbourhood_features,spatial_confounding
from src.validation.gates import simulation_allowed,evaluate_gates
from src.representation.uncertainty import build_uncertainty
from src.synthetic.mechanisms import clonal_selection,transcriptional_plasticity,selection_plus_plasticity

class ExtendedTests(unittest.TestCase):
 def test_checksum(self):
  with tempfile.NamedTemporaryFile() as f:
   f.write(b'x'); f.flush(); self.assertEqual(sha256(f.name),hashlib.sha256(b'x').hexdigest())
 def test_distribution(self):
  d=TumourStateDistribution.from_states(np.array([[0,1],[2,3]]),['a','b']); self.assertEqual(d.mean,[1,2]); self.assertGreater(d.entropy,0)
 def test_spatial(self):
  g=radius_graph([[0,0],[0,1],[10,10]],2); self.assertEqual(g[0],(1,)); self.assertEqual(neighbourhood_context(g,['a','b','c'],2).niche,'UNRESOLVED')
  def test_spatial_missing_coordinates(self):
   with self.assertRaises(ValueError): validate_coordinates([[0,0],[float('nan'),1]])
 def test_observation(self):
  m=LinearGaussianObservation('bulk',np.eye(2),np.ones(2)); self.assertEqual(m.reconstruction_error([1,2],[1,2]),0)
 def test_posterior(self):
  mean,var=gaussian_posterior([0],[1],[2],[1]); self.assertAlmostEqual(float(mean[0]),1); self.assertAlmostEqual(float(var[0]),.5)
 def test_clone_assignment(self): self.assertEqual(assign_clone({'a':.1,'b':.5},.2),'a')
 def test_threshold_sensitivity(self):
  import pandas as pd
  self.assertEqual(cell_threshold_sensitivity(pd.DataFrame({'patient_id':['a','b'],'n_cells':[20,60]}),[50])[50],['b'])
 def test_negative_controls(self): self.assertEqual(len(negative_controls()),6)
 def test_permutation_validation(self):
  with self.assertRaises(ValueError): permutation_test([1,2],[0,1],lambda x,y:0,10)
 def test_adjustment_reports_three_results(self):
  r=three_way_adjustment([1,2,3,4],[1,2,3,4],[0,1,0,1]); self.assertTrue({'unadjusted','adjusted','residualized'} <= r.keys())
 def test_perturbations_cannot_claim_validation(self):
  with self.assertRaises(ValueError): Perturbation('pathway inhibition','IFN','down','VALIDATED')
 def test_synthetic_known_proportions(self):
  d=clonal_selection(); self.assertAlmostEqual(np.mean(d['baseline_clone']==0),.7); self.assertAlmostEqual(np.mean(d['followup_clone']==0),.2)
  p=transcriptional_plasticity(); self.assertTrue(np.array_equal(p['baseline_clone'],p['followup_clone']))
  self.assertEqual(selection_plus_plasticity()['truth'],'selection_plus_plasticity')
 def test_uncertainty_components(self):
  a=build_uncertainty(1,10,.5,1); b=build_uncertainty(1,500,1,1)
  self.assertGreater(a.sampling_uncertainty,b.sampling_uncertainty); self.assertGreater(a.coverage_penalty,b.coverage_penalty)
  self.assertEqual(build_uncertainty(1,500,1,0,('CNV',)).overall_quality,'INSUFFICIENT')
 def test_probabilistic_transition_distribution(self):
  rng=np.random.default_rng(2); z=rng.normal(size=(12,2)); t=np.array(['A']*12); y=z+.5+rng.normal(scale=.1,size=(12,2)); dt=np.arange(1,13); pid=np.array([f'p{i}' for i in range(12)])
  m=__import__('src.dynamics.discrete',fromlist=['ProbabilisticTransition']).ProbabilisticTransition().fit(z,y,t,dt,pid)
  d=m.predict_distribution(z[0],'A',1); self.assertEqual(d['mean'].shape,(1,2)); self.assertEqual(m.sample_next(z[0],'A',1,3).shape,(3,2)); self.assertTrue(np.isfinite(m.log_likelihood(z[0],y[0],'A',1)))
  with self.assertRaises(ValueError): m.predict_mean(z[0],'UNKNOWN',1)
 def test_probabilistic_duplicate_patient_rejected(self):
  m=ProbabilisticTransition()
  with self.assertRaises(ValueError): m.fit(np.zeros((5,1)),np.ones((5,1)),['A']*5,[1]*5,['p','p','q','r','s'])
 def test_spatial_summary_and_confounding(self):
  f=neighbourhood_features(np.array([[0.],[1.],[3.]]),[(1,),(0,2),()],["a","a","b"]); self.assertTrue(f[2]['isolated'])
  self.assertIn('identifiability_status',spatial_confounding([0],[1],{'r':[[0]]},{'r':[[1]]}))
 def test_simulation_is_fail_closed(self):
  self.assertFalse(simulation_allowed({'G1':'PASSED','G2':'PASSED','G9':'PASSED'})); self.assertIn('G10',evaluate_gates({})['gates'])
if __name__=='__main__': unittest.main()
