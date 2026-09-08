import hashlib,tempfile,unittest
from pathlib import Path
import numpy as np
from src.clones.clone_assignment import assign_clone
from src.dynamics.inference import gaussian_posterior
from src.experiments import sha256
from src.observation.model import LinearGaussianObservation
from src.perturbations import Perturbation
from src.representation.latent_state import TumourStateDistribution
from src.spatial.representation import radius_graph,neighbourhood_context
from src.validation.adjustment import three_way_adjustment
from src.validation.patient_level import cell_threshold_sensitivity
from src.validation.permutation import negative_controls

class ExtendedTests(unittest.TestCase):
 def test_checksum(self):
  with tempfile.NamedTemporaryFile() as f:
   f.write(b'x'); f.flush(); self.assertEqual(sha256(f.name),hashlib.sha256(b'x').hexdigest())
 def test_distribution(self):
  d=TumourStateDistribution.from_states(np.array([[0,1],[2,3]]),['a','b']); self.assertEqual(d.mean,[1,2]); self.assertGreater(d.entropy,0)
 def test_spatial(self):
  g=radius_graph([[0,0],[0,1],[10,10]],2); self.assertEqual(g[0],(1,)); self.assertEqual(neighbourhood_context(g,['a','b','c'],2).niche,'UNRESOLVED')
 def test_observation(self):
  m=LinearGaussianObservation('bulk',np.eye(2),np.ones(2)); self.assertEqual(m.reconstruction_error([1,2],[1,2]),0)
 def test_posterior(self):
  mean,var=gaussian_posterior([0],[1],[2],[1]); self.assertAlmostEqual(float(mean[0]),1); self.assertAlmostEqual(float(var[0]),.5)
 def test_clone_assignment(self): self.assertEqual(assign_clone({'a':.1,'b':.5},.2),'a')
 def test_threshold_sensitivity(self):
  import pandas as pd
  self.assertEqual(cell_threshold_sensitivity(pd.DataFrame({'patient_id':['a','b'],'n_cells':[20,60]}),[50])[50],['b'])
 def test_negative_controls(self): self.assertEqual(len(negative_controls()),6)
 def test_adjustment_reports_three_results(self):
  r=three_way_adjustment([1,2,3,4],[1,2,3,4],[0,1,0,1]); self.assertTrue({'unadjusted','adjusted','residualized'} <= r.keys())
 def test_perturbations_cannot_claim_validation(self):
  with self.assertRaises(ValueError): Perturbation('pathway inhibition','IFN','down','VALIDATED')
if __name__=='__main__': unittest.main()
