"""Outcome-independent malignant observability qualification."""
LEVELS=('MALIGNANT_GROUND_TRUTH','MALIGNANT_REPRODUCIBLE_INFERENCE','MALIGNANT_PARTIAL','MALIGNANT_NOT_OBSERVABLE')
def qualify_observability(*, cell_metadata=False, validated_ground_truth=False,
                          cnv_reproducible=False, cnv_pipeline_supported=False,
                          author_annotation=False, author_annotation_type='unknown',
                          author_annotation_provenance='UNKNOWN',
                          author_annotation_validated=False, outcome_used=False,
                          treatment_used=False, epithelial_only=False):
    if outcome_used or treatment_used: return 'MALIGNANT_NOT_OBSERVABLE'
    if validated_ground_truth and cell_metadata: return 'MALIGNANT_GROUND_TRUTH'
    explicit_malignant = (author_annotation and author_annotation_type == 'malignant' and
                          author_annotation_validated and author_annotation_provenance in ('E4','E5'))
    if cell_metadata and explicit_malignant: return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if cell_metadata and cnv_reproducible and cnv_pipeline_supported: return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if cnv_reproducible or author_annotation or epithelial_only or author_annotation_type in ('epithelial','tumour_like','malignant'):
        return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
def phase5_eligible(status): return status in ('MALIGNANT_GROUND_TRUTH','MALIGNANT_REPRODUCIBLE_INFERENCE')
