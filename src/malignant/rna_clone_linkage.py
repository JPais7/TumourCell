"""Fail-closed qualification of public RNA-cell to clone linkage."""
def qualify_rna_clone_linkage(*, same_cell=False, public_inputs=False, method_documented=False,
    malignant_validated=False, clone_assigned=False, direct_ground_truth=False,
    controlled_metadata=False, outcome_used=False, treatment_used=False,
    second_audit_passed=False):
    if outcome_used or treatment_used or controlled_metadata:
        return 'MALIGNANT_NOT_OBSERVABLE'
    if same_cell and public_inputs and malignant_validated and clone_assigned and direct_ground_truth and second_audit_passed:
        return 'MALIGNANT_GROUND_TRUTH'
    if same_cell and public_inputs and method_documented and malignant_validated and clone_assigned and second_audit_passed:
        return 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if same_cell or clone_assigned or malignant_validated:
        return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'

def clone_observability(*, cell_mapping=False, public_inputs=False, method_documented=False, leakage=False):
    if cell_mapping and public_inputs and method_documented and not leakage: return 'VERIFIED'
    if cell_mapping or method_documented: return 'PARTIAL'
    return 'NOT_AVAILABLE'
