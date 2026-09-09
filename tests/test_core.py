import json, tempfile, unittest
from pathlib import Path
import numpy as np
import pandas as pd
from src.representation.encoder import VersionedEncoder
from src.representation.uncertainty import coverage_adjusted_uncertainty
from src.validation.patient_level import patient_bootstrap,patient_pseudobulk,leave_one_patient_out
from src.validation.permutation import permutation_test
from src.validation.leakage import assert_no_holdout_leakage
from src.validation.gates import evaluate_gates
from src.clones.clone_state import identify_mechanisms,NOT_IDENTIFIABLE
from src.dynamics.discrete import transition_matrix
from src.dynamics.transition import GaussianTransition
from src.synthetic.mechanisms import generate,recover

ROOT=Path(__file__).parents[1]
class CoreTests(unittest.TestCase):
 def test_encoder_deterministic_and_missing(self):
  e=VersionedEncoder.from_yaml(ROOT/'configs/representation/latent_state_v1.yaml')
  a=e.encode({'MKI67':2,'TOP2A':1},sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=10,available_modalities=['scRNA-seq'])
  b=e.encode({'MKI67':2,'TOP2A':1},sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=10,available_modalities=['scRNA-seq'])
  self.assertEqual(e.encoder_hash,b.encoder_hash); self.assertEqual(a.features['proliferation'],b.features['proliferation'])
  self.assertIn('CNV',a.missing_modalities); self.assertEqual(a.missing_required_modalities,()); self.assertEqual(a.features['proliferation'].qc,'INSUFFICIENT')
 def test_uncertainty_increases(self):
  self.assertGreater(coverage_adjusted_uncertainty(1,10,.5,.5),coverage_adjusted_uncertainty(1,500,1,1))
 def test_patient_aggregation(self):
  f=pd.DataFrame({'cell_id':['a','b'],'sample_id':['s','s'],'patient_id':['p','p'],'cohort_id':['c','c'],'x':[1,2]})
  self.assertEqual(patient_pseudobulk(f,['x']).x.iloc[0],3)
 def test_bootstrap_seed_and_lopo(self):
  self.assertEqual(patient_bootstrap([1,2,3],100,7),patient_bootstrap([1,2,3],100,7)); self.assertEqual(len(leave_one_patient_out(['a','b'],len)),2)
 def test_permutation_seed(self):
  stat=lambda x,y:np.mean(x[y==1])-np.mean(x[y==0])
  self.assertEqual(permutation_test(np.arange(6),[0,0,0,1,1,1],stat,100,3),permutation_test(np.arange(6),[0,0,0,1,1,1],stat,100,3))
 def test_no_clone_is_not_identifiable(self):
  x=identify_mechanisms(has_clone_data=True,assignment_quality=.95,clones_by_timepoint={'Pre':{'A':.7,'B':.3},'On':{'A':.2,'B':.8}},temporal_overlap=True,sampling_qc=True); self.assertEqual(x['selection_status'],'ESTIMABLE')
  y=identify_mechanisms(has_clone_data=True); self.assertEqual(y['selection_status'],NOT_IDENTIFIABLE); self.assertEqual(y['plasticity_status'],NOT_IDENTIFIABLE)
 def test_transition_and_treatment(self):
  self.assertTrue(np.allclose(transition_matrix([0,0],[1,1],2).sum(1),1))
  m=GaussianTransition().fit([[0],[0]],[[1],[2]],['a','b']); self.assertEqual(float(m.predict([0],'a')[0]),1)
 def test_gate_non_cascade(self):
  self.assertFalse(evaluate_gates({'G1':'PASSED','G2':'PASSED'})['simulation_enabled'])
 def test_leakage(self):
  with self.assertRaises(ValueError): assert_no_holdout_leakage([{'role':'holdout','use':'encoder_fitting'}])
 def test_synthetic_truths(self):
  for mechanism in ('selection','plasticity','selection_plus_plasticity','sampling_depth','spatial_environment'):
   self.assertEqual(recover(generate(mechanism,1000,4)),mechanism)
 def test_config_schema_fields(self):
  e=VersionedEncoder.from_yaml(ROOT/'configs/representation/latent_state_v1.yaml'); d=e.as_dict()
  self.assertIn('encoder_hash',d); self.assertEqual(len(d['programs']),15)
 def test_raw_counts_transform_and_modalities(self):
  e=VersionedEncoder.from_yaml(ROOT/'configs/representation/latent_state_v1.yaml')
  raw={'MKI67':10,'TOP2A':20,'UBE2C':30}
  a=e.encode_counts(raw,sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=100,available_modalities=['scRNA-seq'])
  total=sum(raw.values()); transformed={k:np.log1p(v/total*1e6) for k,v in raw.items()}
  b=e.encode(transformed,sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=100,available_modalities=['scRNA-seq'])
  self.assertAlmostEqual(a.features['proliferation'].value,b.features['proliferation'].value)
  self.assertEqual(a.features['proliferation'].modality_coverage,1.0)
 def test_incompatible_modality_fails(self):
  e=VersionedEncoder.from_yaml(ROOT/'configs/representation/latent_state_v1.yaml')
  with self.assertRaises(ValueError): e.encode({'MKI67':1},sample_id='s',patient_id='p',cohort_id='c',timepoint='t',n_cells=50,available_modalities=['foobar'])
if __name__=='__main__': unittest.main()
