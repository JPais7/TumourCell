"""Population and cohort priors."""
import numpy as np
def population_prior(states):
    x=np.asarray(states,float)
    return {"mean":x.mean(0),"variance":x.var(0,ddof=1),"n_patients":len(x)}
def cohort_conditioned_prior(states,cohorts,target):
    mask=np.asarray(cohorts)==target
    if not mask.any(): raise ValueError("target cohort absent")
    return population_prior(np.asarray(states)[mask])
