"""Population-prior → cohort-conditioned → patient-informed posterior."""
import numpy as np
def gaussian_posterior(prior_mean,prior_variance,observation,observation_variance):
    pm,pv,o,ov=map(lambda x:np.asarray(x,float),(prior_mean,prior_variance,observation,observation_variance))
    precision=1/pv+1/ov
    return (pm/pv+o/ov)/precision,1/precision
