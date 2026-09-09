"""Outcome-independent malignant observability qualification."""
LEVELS=('MALIGNANT_GROUND_TRUTH','MALIGNANT_REPRODUCIBLE_INFERENCE','MALIGNANT_PARTIAL','MALIGNANT_NOT_OBSERVABLE')
def qualify_observability(*, cell_metadata=False, validated_ground_truth=False, cnv_reproducible=False, author_annotation=False, outcome_used=False, epithelial_only=False):
    if outcome_used: return 'MALIGNANT_NOT_OBSERVABLE'
    if validated_ground_truth and cell_metadata: return 'MALIGNANT_GROUND_TRUTH'
    if cell_metadata and (cnv_reproducible or author_annotation): return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if cnv_reproducible or author_annotation or epithelial_only: return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
def phase5_eligible(status): return status in ('MALIGNANT_GROUND_TRUTH','MALIGNANT_REPRODUCIBLE_INFERENCE')
