# Phase 4.1 scientific-integrity audit

This pass hardens the existing Phase 4 implementation without changing Phase 1–3 definitions, hashes or historical outputs.

Fixed:

- Leave-one-patient-out transition validation now fits each fold independently and uses that fold's covariance for marginal predictive intervals. Fold predictions include patient, treatment, `delta_t`, mean, covariance, residual and fold ID.
- Directional accuracy compares `z_next - z_current`, per dimension and overall, rather than signs of absolute states.
- The legacy boolean mechanism classifier is now a deprecated fail-closed wrapper around `ObservedShiftAssessment`. Clone evidence alone cannot support selection; expression change alone cannot support plasticity.
- Spatial decomposition now reports global shift, within-region shift and descriptive between-region composition shift when regions overlap, otherwise `NOT_IDENTIFIABLE` with reasons.
- A single `SIMULATION_PREREQUISITE_GATES` policy controls gate evaluation. `simulation_allowed()` remains permanently fail-closed in Phase 4 and clinical recommendations remain disabled.

Status distinction:

- Probabilistic dynamics: IMPLEMENTED, UNIT_TESTED, SYNTHETICALLY_VALIDATED; real-cohort validation NOT_IMPLEMENTED.
- Mechanism decomposition: IMPLEMENTED, UNIT_TESTED, SYNTHETICALLY_VALIDATED; real-cohort identifiability NOT_IDENTIFIABLE.
- Spatial observation: IMPLEMENTED, UNIT_TESTED; causal spatial validation NOT_IDENTIFIABLE/NOT_TESTED.
- Observation model: INTERFACE_ONLY; not biologically validated.

`pytest -q` passes 32 tests. The project still obeys: association is not causality; prediction is not dynamics; dynamics is not counterfactual validity; synthetic validation is not biological validation.
