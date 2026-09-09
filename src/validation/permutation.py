"""Patient-level randomization controls."""
from __future__ import annotations
import numpy as np

def permutation_test(values, labels, statistic, iterations=10000, seed=0):
    x=np.asarray(values); y=np.asarray(labels); rng=np.random.default_rng(seed)
    if iterations < 100: raise ValueError("at least 100 permutations required")
    if len(x) != len(y) or len(x) < 4: raise ValueError("values and labels require >=4 aligned observations")
    if len(np.unique(y)) != 2: raise ValueError("labels must contain exactly two groups")
    observed=float(statistic(x,y)); null=np.empty(iterations)
    for i in range(iterations): null[i]=statistic(x,rng.permutation(y))
    return {"observed": observed, "effect_size": observed, "p_two_sided": float((1+np.sum(np.abs(null)>=abs(observed)))/(iterations+1)), "n_permutations": iterations, "seed": seed, "null_distribution": null.tolist()}

def negative_controls():
    return ("shuffled_treatment","shuffled_response","holdout_cohort","platform_holdout","negative_control_cohort","technical_perturbation")
