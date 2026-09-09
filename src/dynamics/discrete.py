"""Discrete transition counting with smoothing."""
import numpy as np
from dataclasses import dataclass
def transition_matrix(current,next_state,n_states,alpha=1.0):
    counts=np.full((n_states,n_states),alpha,float)
    for a,b in zip(current,next_state): counts[int(a),int(b)]+=1
    return counts/counts.sum(1,keepdims=True)

@dataclass
class ContinuousTransition:
    """Interpretable p(z_next | z_current, treatment), patient-level rows only."""
    coefficients: dict[str, np.ndarray] | None = None
    residual_variance: dict[str, np.ndarray] | None = None
    def fit(self, current_state, next_state, treatment, patient_ids=None):
        z=np.asarray(current_state,float); y=np.asarray(next_state,float); t=np.asarray(treatment)
        if z.ndim!=2 or y.shape!=z.shape or len(t)!=len(z): raise ValueError("unaligned continuous transitions")
        if patient_ids is not None and len(np.unique(patient_ids)) < len(z):
            raise ValueError("duplicate patient rows: aggregate to patient before fitting")
        self.coefficients={}; self.residual_variance={}
        for label in np.unique(t):
            mask=t==label
            if mask.sum()<2: raise ValueError(f"at least two patient transitions required for treatment {label}")
            X=np.column_stack([np.ones(mask.sum()),z[mask]])
            beta=np.linalg.lstsq(X,y[mask],rcond=None)[0]; self.coefficients[str(label)]=beta
            self.residual_variance[str(label)]=np.var(y[mask]-X@beta,axis=0,ddof=1)
        return self
    def predict(self,current_state,treatment):
        if self.coefficients is None or str(treatment) not in self.coefficients: raise ValueError("unfitted or unseen treatment")
        z=np.asarray(current_state,float); X=np.column_stack([np.ones(len(z)),z]) if z.ndim==2 else np.r_[1,z][None,:]
        return X@self.coefficients[str(treatment)]
