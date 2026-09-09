def qualify_literature_clone(*, direct_mapping=False, public_inputs=False, genomic_evidence=False, malignant_validated=False, mapping_reproducible=False, second_audit=False, controlled=False, ambiguous=False, outcome_used=False, treatment_used=False):
    if controlled or outcome_used or treatment_used or ambiguous: return 'MALIGNANT_NOT_OBSERVABLE'
    if direct_mapping and public_inputs and genomic_evidence and malignant_validated and second_audit:
        return 'MALIGNANT_GROUND_TRUTH' if mapping_reproducible else 'MALIGNANT_REPRODUCIBLE_INFERENCE'
    if direct_mapping or genomic_evidence: return 'MALIGNANT_PARTIAL'
    return 'MALIGNANT_NOT_OBSERVABLE'
def reconstruct_mapping(*, n_rna_cells=0, n_cells_with_clone=0, unique=True, public_inputs=True, reference=True):
    if not (n_rna_cells and unique and public_inputs and reference): return {'successful': False, 'fraction': 0.0}
    return {'successful': n_cells_with_clone > 0, 'fraction': n_cells_with_clone / n_rna_cells}
