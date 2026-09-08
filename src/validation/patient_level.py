"""Central cell→sample→patient→cohort aggregation and resampling."""
from __future__ import annotations
import numpy as np
import pandas as pd

REQUIRED_HIERARCHY = ("cell_id", "sample_id", "patient_id", "cohort_id")

def validate_hierarchy(frame: pd.DataFrame) -> None:
    missing = set(REQUIRED_HIERARCHY) - set(frame.columns)
    if missing: raise ValueError(f"missing hierarchy columns: {sorted(missing)}")
    for child, parent in zip(REQUIRED_HIERARCHY, REQUIRED_HIERARCHY[1:]):
        if (frame.groupby(child)[parent].nunique() > 1).any(): raise ValueError(f"{child} maps to multiple {parent}s")

def patient_pseudobulk(frame: pd.DataFrame, value_columns: list[str]) -> pd.DataFrame:
    validate_hierarchy(frame)
    samples = frame.groupby(["cohort_id","patient_id","sample_id"], as_index=False)[value_columns].sum()
    return samples.groupby(["cohort_id","patient_id"], as_index=False)[value_columns].mean()

def patient_bootstrap(values, iterations=10000, seed=0, confidence=.95):
    x=np.asarray(values,float); rng=np.random.default_rng(seed)
    means=np.asarray([rng.choice(x,len(x),replace=True).mean() for _ in range(iterations)])
    a=(1-confidence)/2
    return float(x.mean()), tuple(float(v) for v in np.quantile(means,[a,1-a]))

def leave_one_patient_out(patient_ids, estimator):
    ids=list(dict.fromkeys(patient_ids))
    return {p: estimator([x for x in ids if x != p]) for p in ids}

def cell_threshold_sensitivity(frame, thresholds, patient="patient_id", cells="n_cells"):
    return {int(t): sorted(frame.loc[frame[cells] >= t, patient].unique().tolist()) for t in thresholds}
