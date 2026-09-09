# Phase 4.3.2 — Selection, threshold and model-selection leakage audit

Repository inspection found no learned feature-selection, threshold-learning, model-family selection, hyperparameter tuning or nested temporal model-selection API in the current Phase 4 forecasting path. These statuses are therefore `NOT_IMPLEMENTED`, not `UNIT_TESTED`. The historical Phase 1 NMF/gene-discovery scripts are frozen outcome-blind artefacts and were not altered or retuned here.

The implemented adversarial tests use an extreme future-only feature and demonstrate that a historical-only scaler does not use it, while deliberately fitting on historical plus future data changes the parameters. This tests the membership contract of an existing preprocessing component; it does not manufacture a feature selector or model selector.

Patient-level LOO, strict temporal history/target validation and frozen-encoder immutability remain separately tested. Joint patient-plus-temporal generalization remains `NOT_IMPLEMENTED`. Real-cohort forecasting remains `NOT_IMPLEMENTED`, Phase 5 remains blocked, and simulation/clinical recommendations remain disabled.

Passing these software tests does not establish biological validity, treatment causality, real-tumour forecasting performance or clinical utility.
