"""Cell-to-genome linkage qualification; fail closed without exact mapping."""
def qualify_linkage(*, same_cell=False, malignant_validated=False, public_inputs=False, method_documented=False, reference_available=False, direct_ground_truth=False, outcome_used=False, treatment_used=False):
    if outcome_used or treatment_used: return 'MALIGNANT_NOT_OBSERVABLE'
    if same_cell and malignant_validated and public_inputs and direct_ground_truth: return 'MALIGNANT_GROUND_TRUTH'
    if same_cell and public_inputs and method_documented and reference_available and malignant_validated: return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if same_cell or malignant_validated: return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
def phase5_eligible_linkage(status): return status in ('MALIGNANT_GROUND_TRUTH','MALIGNANT_REPRODUCIBLE_INFERENCE')
