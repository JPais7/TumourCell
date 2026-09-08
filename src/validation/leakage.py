"""Fail-closed data-role checks for frozen validation."""
FORBIDDEN_HOLDOUT_USES={"feature_selection","programme_discovery","encoder_fitting","hyperparameter_tuning","threshold_selection","model_selection"}
def assert_no_holdout_leakage(records):
    bad=[r for r in records if r.get("role") in {"holdout","external_validation"} and r.get("use") in FORBIDDEN_HOLDOUT_USES]
    if bad: raise ValueError(f"holdout leakage detected: {bad}")
    return True
