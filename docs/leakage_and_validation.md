# Leakage and validation contract

Frozen program definitions, encoder fitting, treatment labels, response labels, threshold selection and model selection use discovery data only. Holdout patients and validation cohorts are excluded from these operations and cannot alter state names. Longitudinal dynamics use patient-level transition rows; random cell-level splits are prohibited. Spatial parameters are fitted only inside the training region/cohort. Clinical outcomes are never used to define P8 or State_01–State_04. External prediction does not imply dynamics, causality or counterfactual validity.

Synthetic truth is never used to tune thresholds on real cohorts. Frozen encoders remain frozen during external validation. A spatial validation region cannot influence training. Patient IDs are checked for disjointness across temporal folds. These are integrity controls, not evidence that a biological mechanism has been validated.
