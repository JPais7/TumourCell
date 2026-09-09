# Phase 4.4 — Real-cohort temporal forecasting benchmark

The frozen GSE246613 pseudobulk table was used without redefining P8 or State_01–State_04. It contains 99 pseudobulks from 34 patients. After patient/time ordering, 65 patient-level transitions were available; response labels R1/R2/NR were not predictors.

In patient-held-out forecasting, persistence was better than the probabilistic transition model: persistence MAE 0.927 and RMSE 1.351 versus model MAE 1.122 and RMSE 1.499. The model directional accuracy was 0.554, with per-feature values [0.585, 0.569, 0.508, 0.585, 0.523]. The 95% predictive interval coverage was 0.846 for the model and 0.88 for the descriptive persistence comparison. These are observational predictive metrics, not evidence of treatment causality or biological mechanism.

Strict forward-temporal validation was not estimable at the first boundary: historical Base transitions do not provide a sufficient PD1-conditioned training model, and pooling with future PD1 transitions would violate the temporal contract. Joint patient-plus-temporal validation is not implemented. The benchmark is therefore a real-cohort software benchmark with limited externality, not a validated tumour dynamics model.

The result is `MODEL_WORSE_ON_MAE_AND_RMSE`; no post-hoc tuning was performed. Phase 5 remains blocked, simulation is disabled, and clinical recommendations remain disabled.
