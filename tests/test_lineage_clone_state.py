from src.dynamics.lineage_clone_state import *
def test_patient_not_clone(): assert qualify_lineage_clone_state(rna_cells=True,repeated_phases=True,state_representation=True)=='INSUFFICIENT'
def test_barcode_mapping_required(): assert not selection_plasticity_testable(clone_frequencies=False,within_clone_states=True,repeated_phases=True)
def test_repeated_phases_required(): assert not selection_plasticity_testable(clone_frequencies=True,within_clone_states=True,repeated_phases=False)
def test_same_cell_lineage_not_required(): assert selection_plasticity_testable(clone_frequencies=True,within_clone_states=True,repeated_phases=True)
def test_valid_clonmapper(): assert qualify_lineage_clone_state(rna_cells=True,barcode_mapping=True,repeated_phases=True,state_representation=True)=='LINEAGE_CLONE_STATE_OBSERVABLE'
def test_genomic_clone_not_assumed(): assert qualify_lineage_clone_state(rna_cells=True,barcode_mapping=True,repeated_phases=True,state_representation=True,genomic_clone=False)=='LINEAGE_CLONE_STATE_OBSERVABLE'
def test_controlled_egas(): assert qualify_lineage_clone_state(rna_cells=True,barcode_mapping=True,repeated_phases=True,state_representation=True,controlled_inputs=True,public_inputs=False)=='CONTROLLED'
def test_treatment_leakage(): assert qualify_lineage_clone_state(treatment_used=True)=='NOT_OBSERVABLE'
def test_outcome_leakage(): assert qualify_lineage_clone_state(outcome_used=True)=='NOT_OBSERVABLE'
def test_phase5_locked(): assert phase5_locked()=={'phase5_ready':False,'simulation_allowed':False,'clinical_recommendations_enabled':False}
