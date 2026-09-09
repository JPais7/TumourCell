from src.malignant.literature_clone import *
def test_patient_not_cell(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_sample_not_cell(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_tumour_not_cell(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_cluster_not_clone(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_publication_only(): assert qualify_literature_clone()=='MALIGNANT_NOT_OBSERVABLE'
def test_clone_table_no_barcode(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_rna_no_genome(): assert qualify_literature_clone(direct_mapping=True)=='MALIGNANT_PARTIAL'
def test_genome_no_mapping(): assert qualify_literature_clone(genomic_evidence=True)=='MALIGNANT_PARTIAL'
def test_similar_barcodes(): assert qualify_literature_clone(direct_mapping=True)=='MALIGNANT_PARTIAL'
def test_row_order(): assert not reconstruct_mapping(n_rna_cells=10,n_cells_with_clone=10,public_inputs=False)['successful']
def test_controlled(): assert qualify_literature_clone(controlled=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_reference(): assert not reconstruct_mapping(n_rna_cells=10,n_cells_with_clone=5,reference=False)['successful']
def test_epithelial(): assert qualify_literature_clone(direct_mapping=True,public_inputs=True,genomic_evidence=True,second_audit=True)!='MALIGNANT_GROUND_TRUTH'
def test_treatment(): assert qualify_literature_clone(treatment_used=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_outcome(): assert qualify_literature_clone(outcome_used=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_future_audit(): assert qualify_literature_clone(second_audit=False)=='MALIGNANT_NOT_OBSERVABLE'
def test_e4(): assert qualify_literature_clone(direct_mapping=True,public_inputs=True,genomic_evidence=True,malignant_validated=True,mapping_reproducible=False,second_audit=True)=='MALIGNANT_REPRODUCIBLE_INFERENCE'
def test_e5(): assert qualify_literature_clone(direct_mapping=True,public_inputs=True,genomic_evidence=True,malignant_validated=True,mapping_reproducible=True,second_audit=True)=='MALIGNANT_GROUND_TRUTH'
def test_partial(): assert qualify_literature_clone(direct_mapping=True)=='MALIGNANT_PARTIAL'
def test_ambiguous(): assert qualify_literature_clone(direct_mapping=True,ambiguous=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_many_to_one(): assert not reconstruct_mapping(n_rna_cells=10,n_cells_with_clone=10,unique=False)['successful']
def test_clone_not_lineage(): assert True
def test_repeated_without_mapping(): assert not reconstruct_mapping(n_rna_cells=10,n_cells_with_clone=0)['successful']
def test_phase5_locked(): assert True
def test_zero_cells(): assert not reconstruct_mapping()['successful']
def test_missing_public(): assert not reconstruct_mapping(n_rna_cells=1,n_cells_with_clone=1,public_inputs=False)['successful']
