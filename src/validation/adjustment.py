"""Triangulate unadjusted, adjusted and residualized estimates."""
import numpy as np
def three_way_adjustment(outcome,signal,covariates):
    y=np.asarray(outcome,float); x=np.asarray(signal,float); c=np.asarray(covariates,float)
    if c.ndim==1:c=c[:,None]
    design=np.column_stack([np.ones(len(y)),x,c]); adjusted=float(np.linalg.lstsq(design,y,rcond=None)[0][1])
    nuisance=np.column_stack([np.ones(len(y)),c]); residual=x-nuisance@np.linalg.lstsq(nuisance,x,rcond=None)[0]
    residualized=float(np.linalg.lstsq(np.column_stack([np.ones(len(y)),residual]),y,rcond=None)[0][1])
    unadjusted=float(np.linalg.lstsq(np.column_stack([np.ones(len(y)),x]),y,rcond=None)[0][1])
    return {"unadjusted":unadjusted,"adjusted":adjusted,"residualized":residualized,
            "adjustment_removed_potential_biological_signal":abs(adjusted)<.5*abs(unadjusted) if unadjusted else False}
