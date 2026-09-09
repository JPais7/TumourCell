"""Fail-closed exhaustive cell-to-genome linkage qualification."""
def qualify_exhaustive_candidate(*, same_cell=False, public_mapping=False, malignant_validated=False, direct_ground_truth=False, cnv_supported=False, reference_available=False, second_audit_passed=False, publication_only=False, outcome_used=False, treatment_used=False):
    if outcome_used or treatment_used: return 'MALIGNANT_NOT_OBSERVABLE'
    if publication_only: return 'MALIGNANT_PARTIAL'
    if not second_audit_passed and (same_cell or public_mapping): return 'MALIGNANT_PARTIAL'
    if same_cell and public_mapping and malignant_validated and direct_ground_truth and second_audit_passed: return 'MALIGNANT_GROUND_TRUTH'
    if same_cell and public_mapping and malignant_validated and cnv_supported and reference_available and second_audit_passed: return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if same_cell or malignant_validated or cnv_supported: return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
