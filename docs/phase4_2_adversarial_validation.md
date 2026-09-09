# Phase 4.2 adversarial scientific validation

Phase 4.2 attempted to break the computational infrastructure before any real longitudinal dynamics are fitted. The tests target optimistic validation, patient leakage, duplicate transitions, treatment imbalance, singular covariance, wrong directional metrics, insufficient evidence for mechanism attribution, spatial decomposition errors, malformed coordinates and gate bypasses.

`pytest -q` collected and passed 45 tests; no tests were skipped or xfailed. Leave-one-patient-out folds expose train/test patient identifiers and use covariance from the fitted training fold. Directional accuracy is computed from `z_next - z_current`. Gaussian prediction intervals remain a model assumption, not a claim of calibrated uncertainty on biological data.

The mechanism assessment distinguishes `SUPPORTED`, `NOT_SUPPORTED`, `INDETERMINATE` and `NOT_IDENTIFIABLE`. Clone evidence alone cannot support selection, expression change alone cannot support plasticity, and spatial correlation is never interpreted as a causal niche. Missing overlap, empty regions and sampling/technical confounding fail closed or remain indeterminate.

The gate policy is exhaustive and fail-closed: `simulation_allowed()` remains `False` even if a synthetic evidence dictionary marks every gate as passed. Clinical recommendations remain disabled. Passing these synthetic adversarial tests demonstrates only that the implementation obeys its specified rules; it does not establish tumour biology, real-cohort dynamics, causality or Digital Twin validity.

Historical Phase 1–3 results, frozen definitions and hashes were not modified. The strongest supported conclusion is that the machinery is adversarially tested while biological validation remains unestablished.
