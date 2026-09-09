"""Conservative ClonMapper reconstruction and selection/plasticity gates."""
def qualify_clone_state_reconstruction(*, rna_cells=0, lineage_mapping=False, mapping_reproducible=False, state_available=False, repeated_phases=False, controlled=False, row_order=False, namespace_match=True, ambiguous=False):
    if controlled: return 'CONTROLLED_DATA_REQUIRED'
    if row_order or ambiguous or not namespace_match: return 'RECONSTRUCTION_FAILED'
    if rna_cells and lineage_mapping and mapping_reproducible and state_available: return 'RECONSTRUCTION_VERIFIED'
    if rna_cells and (lineage_mapping or state_available): return 'RECONSTRUCTION_PARTIAL'
    return 'RECONSTRUCTION_FAILED'
def qualify_selection_plasticity_evidence(*, clone_frequency=False, within_clone_state=False, baselines_evaluated=False, uncertainty=False, leakage=False):
    if leakage or not baselines_evaluated or not uncertainty: return 'INSUFFICIENT_DATA'
    if clone_frequency and within_clone_state: return 'MIXED_EVIDENCE'
    if clone_frequency: return 'SELECTION_SUPPORTED'
    if within_clone_state: return 'WITHIN_CLONE_STATE_REDISTRIBUTION_SUPPORTED'
    return 'INDETERMINATE'
def decomposition_status(*, mapped_cells=0, repeated_phases=False, clone_state_tables=False):
    return 'IDENTIFIABLE' if mapped_cells and repeated_phases and clone_state_tables else 'NOT_IDENTIFIABLE'
