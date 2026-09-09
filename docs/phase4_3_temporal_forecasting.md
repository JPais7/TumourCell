# Phase 4.3 — Temporal forecasting and leakage audit

Phase 4.3 adds a separate temporal contract. A forecast at time `t` may use only rows with time strictly before `t`; future states, outcomes, treatment responses, preprocessing parameters, feature selections, thresholds and covariances are excluded. Patient generalisation, temporal generalisation and joint patient-plus-temporal generalisation remain distinct concepts.

The implementation provides strict temporal splitting, duplicate and ordering checks, transition construction with observed irregular `delta_t`, and a Phase 5 readiness assessment. It does not fit real tumour longitudinal dynamics. The existing patient-level LOO validator remains a patient-generalisation procedure; it has not been relabelled as temporal validation.

Synthetic tests cover stationary linear dynamics, irregular intervals, future-row exclusion, insufficient history, duplicate transitions and malformed timestamps. These tests establish computational behaviour only. Predictive MAE/RMSE, log likelihood, marginal predictive interval coverage and delta-based directional accuracy remain evaluation metrics, not evidence of biological validity or treatment causality.

The current status is `temporal_forecasting = SYNTHETICALLY_VALIDATED`, `real_cohort_forecasting = NOT_IMPLEMENTED`, and `phase5_ready = false`. No synthetic truth tunes real-cohort thresholds. Treatment is an observational conditioning variable only. Simulation, counterfactuals and clinical recommendations remain disabled.
