# Phase 4.3.1 — Temporal leakage hardening

The temporal contract now exposes a `ForecastExample` with explicit history, target, prediction time and observed `delta_t`. It fails closed when history is empty, target time is mismatched, timestamps are duplicated/non-monotonic, or future rows would enter historical fitting.

The hardening tests demonstrate that adding future rows cannot change a historical scaler, a frozen encoder hash, or the covariance used by a model fitted only on historical transitions. These tests inspect fitting membership and parameters, not merely predictive performance. Feature-selection, threshold-learning and nested model-selection APIs do not exist in this repository and are therefore recorded as `NOT_IMPLEMENTED`, not assumed safe.

Patient-level LOO remains a patient-generalization procedure. The new history→target contract is temporal generalization. Joint patient-plus-temporal validation is `NOT_IMPLEMENTED`. No real tumour cohort was fitted. Phase 5 remains blocked, `simulation_allowed()` remains false, and clinical recommendations remain disabled.

The distinction remains explicit: a forecast can be computationally leakage-controlled without establishing biological mechanism, treatment causality or counterfactual validity.
