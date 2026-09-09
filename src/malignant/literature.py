"""Evidence checks for literature-linked cell-to-genome candidates."""
def qualify_literature_candidate(*, same_cell=False, public_mapping=False, malignant_validated=False, cnv_supported=False, reference_available=False, direct_ground_truth=False, outcome_used=False, treatment_used=False):
    if outcome_used or treatment_used: return 'MALIGNANT_NOT_OBSERVABLE'
    if same_cell and public_mapping and malignant_validated and direct_ground_truth: return 'MALIGNANT_GROUND_TRUTH'
    if same_cell and public_mapping and malignant_validated and cnv_supported and reference_available: return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if same_cell or malignant_validated or cnv_supported: return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
