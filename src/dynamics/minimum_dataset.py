"""Fail-closed gates for clone/state longitudinal mechanistic observability."""

def qualify_mechanistic_observability(*, malignant_identity=False, clone_mapping=False,
    longitudinal=False, state_representation=False, treatment_annotated=False,
    outcome_used=False, treatment_used=False, future_baseline=False,
    same_cell_lineage=False, independent_validation=False, persistence_beaten=False,
    observation_uncertainty=False):
    if outcome_used or treatment_used or future_baseline:
        return False
    return all((malignant_identity, longitudinal, state_representation))

def selection_identifiable(*, malignant_identity=False, clone_mapping=False,
    repeated_samples=False, public_mapping=False, leakage=False):
    return bool(malignant_identity and clone_mapping and repeated_samples and public_mapping and not leakage)

def plasticity_identifiable(*, malignant_identity=False, clone_mapping=False,
    repeated_samples=False, within_clone_state=False, same_cell_lineage=False, leakage=False):
    # Same-cell lineage is not required for population-level within-clone comparison.
    return bool(malignant_identity and clone_mapping and repeated_samples and within_clone_state and not leakage)

def phase5a_eligible(*, malignant_e4_or_stronger=False, g1_transport=False,
    g2_external_prediction=False, g3_treatment_transition=False,
    partial_clone_state=False, leakage=False, persistence_beaten=False,
    independent_validation=False, observation_uncertainty=False):
    return all((malignant_e4_or_stronger, g1_transport, g2_external_prediction,
                g3_treatment_transition, partial_clone_state, not leakage,
                persistence_beaten, independent_validation, observation_uncertainty))
