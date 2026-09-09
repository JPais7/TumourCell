from pathlib import Path
from src.validation.real_cohort import frozen_scores,forward_temporal_readiness
ROOT=Path(__file__).parents[1]
def test_gse246613_frozen_representation_has_coverage():
 x,z,names,c=frozen_scores(ROOT/'data/processed/GSE246613/malignant_pseudobulk_counts.npz',ROOT/'results/phase3/state_atlas_v1.json',ROOT/'results/phase2/p8_definition_frozen.json')
 assert len(set(x['patient']))==34 and z.shape[0]==99 and min(c.values())==1.0
def test_forward_temporal_boundary_is_conservative():
 r=forward_temporal_readiness(ROOT/'data/processed/GSE246613/malignant_pseudobulk_counts.npz')
 assert r['status']=='NOT_ESTIMABLE'
