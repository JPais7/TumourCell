"""Fail-closed qualification for experimental lineage-barcode clone-state data."""
def qualify_lineage_clone_state(*, rna_cells=False, barcode_mapping=False, repeated_phases=False,
    state_representation=False, genomic_clone=False, public_inputs=False,
    controlled_inputs=False, outcome_used=False, treatment_used=False):
    if outcome_used or treatment_used: return 'NOT_OBSERVABLE'
    if not (rna_cells and barcode_mapping and repeated_phases and state_representation): return 'INSUFFICIENT'
    if controlled_inputs and not public_inputs: return 'CONTROLLED'
    return 'LINEAGE_CLONE_STATE_OBSERVABLE'

def selection_plasticity_testable(*, clone_frequencies=False, within_clone_states=False,
    repeated_phases=False, same_cell_lineage=False):
    return bool(clone_frequencies and within_clone_states and repeated_phases)

def phase5_locked():
    return {'phase5_ready': False, 'simulation_allowed': False, 'clinical_recommendations_enabled': False}
