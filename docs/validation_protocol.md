# Validation protocol

All inference and resampling are patient-level. Required controls are shuffled treatment, shuffled response, cohort holdout, platform holdout, negative-control cohort and technical perturbation. Leave-one-patient-out and cell-threshold sensitivity are mandatory. Holdouts cannot participate in feature/program discovery, encoder fitting, hyperparameter or threshold selection, or model selection. The historical GSE246613 threshold remains 50 cells.
