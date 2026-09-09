from src.malignant.observability import qualify_observability, phase5_eligible
def test_publication_only_is_not_ground_truth(): assert qualify_observability(author_annotation=True,author_annotation_type='malignant')=='MALIGNANT_PARTIAL'
def test_sample_tumour_label_is_not_cell_malignancy(): assert qualify_observability()== 'MALIGNANT_NOT_OBSERVABLE'
def test_unsupported_cnv_is_partial(): assert qualify_observability(cnv_reproducible=True,cnv_pipeline_supported=False)=='MALIGNANT_PARTIAL'
def test_supported_cnv_is_inference_not_ground_truth(): assert qualify_observability(cell_metadata=True,cnv_reproducible=True,cnv_pipeline_supported=True)=='MALIGNANT_REPRODUCIBLE_INFERENCE'
def test_outcome_and_treatment_fail_closed():
    assert qualify_observability(cell_metadata=True,cnv_reproducible=True,cnv_pipeline_supported=True,outcome_used=True)=='MALIGNANT_NOT_OBSERVABLE'
    assert qualify_observability(cell_metadata=True,cnv_reproducible=True,cnv_pipeline_supported=True,treatment_used=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_partial_never_phase5(): assert phase5_eligible('MALIGNANT_PARTIAL') is False
def test_unknown_semantics_remain_partial_or_unobservable(): assert qualify_observability(author_annotation=True,author_annotation_type='unknown')=='MALIGNANT_PARTIAL'
def test_ground_truth_requires_cell_metadata(): assert qualify_observability(validated_ground_truth=True,cell_metadata=False)=='MALIGNANT_NOT_OBSERVABLE'
