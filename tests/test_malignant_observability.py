from src.malignant.observability import qualify_observability, phase5_eligible
def test_epithelial_only_is_partial(): assert qualify_observability(epithelial_only=True)=='MALIGNANT_PARTIAL'
def test_tumour_sample_is_not_cell_ground_truth(): assert qualify_observability()=='MALIGNANT_NOT_OBSERVABLE'
def test_response_cannot_promote(): assert qualify_observability(cnv_reproducible=True,outcome_used=True)=='MALIGNANT_NOT_OBSERVABLE'
def test_cnv_is_not_ground_truth(): assert qualify_observability(cnv_reproducible=True)=='MALIGNANT_PARTIAL'
def test_downloadable_annotation_can_be_reproducible(): assert qualify_observability(cell_metadata=True,author_annotation=True)=='MALIGNANT_REPRODUCIBLE_INFERENCE'
def test_partial_does_not_unlock(): assert phase5_eligible('MALIGNANT_PARTIAL') is False
