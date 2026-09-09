"""Discrete transition counting with smoothing."""
import numpy as np
from dataclasses import dataclass
import hashlib, json
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

@dataclass
class ProbabilisticTransition:
    """Gaussian p(z_next | z_current, treatment, delta_t), patient-level only."""
    version: str = "probabilistic_transition_v1"
    coefficients: dict | None = None
    residual_covariance: dict | None = None
    feature_names: tuple[str,...] = ()
    treatment_levels: tuple[str,...] = ()
    training_metadata: dict | None = None
    model_hash: str | None = None

    def fit(self, current_state, next_state, treatment, delta_t, patient_id, feature_names=None):
        z=np.asarray(current_state,float); y=np.asarray(next_state,float); t=np.asarray(treatment); dt=np.asarray(delta_t,float); pid=np.asarray(patient_id)
        if z.ndim!=2 or y.shape!=z.shape or len(t)!=len(z) or len(dt)!=len(z) or len(pid)!=len(z): raise ValueError("transition inputs are not aligned")
        if not np.isfinite(z).all() or not np.isfinite(y).all() or not np.isfinite(dt).all(): raise ValueError("NaN/non-finite transition input")
        keys=list(zip(pid,t,dt))
        if len(set(keys)) != len(keys): raise ValueError("duplicate patient × transition rows; aggregate before fitting")
        if np.any(dt<=0): raise ValueError("delta_t must be positive")
        if len(z)<max(5,z.shape[1]+3): raise ValueError("insufficient patient transitions")
        names=tuple(feature_names or [f"z{i}" for i in range(z.shape[1])]); self.feature_names=names
        self.treatment_levels=tuple(sorted(map(str,set(t))))
        self.coefficients={}; self.residual_covariance={}
        for label in self.treatment_levels:
            mask=np.asarray([str(x)==label for x in t]);
            if mask.sum()<3: raise ValueError(f"at least three patient transitions required for treatment {label}")
            # intercept, current z, treatment-specific time slope
            X=np.column_stack([np.ones(mask.sum()),z[mask],dt[mask]])
            beta=np.linalg.lstsq(X,y[mask],rcond=None)[0]; resid=y[mask]-X@beta
            cov=np.atleast_2d(np.cov(resid,rowvar=False)); cov += np.eye(cov.shape[0])*1e-8
            if np.linalg.eigvalsh(cov).min()<=0: raise ValueError("residual covariance is not positive definite")
            self.coefficients[label]=beta.tolist(); self.residual_covariance[label]=cov.tolist()
        self.training_metadata={"n_transitions":len(z),"n_patients":len(set(pid)),"delta_t_min":float(dt.min()),"delta_t_max":float(dt.max()),"unit":"patient"}
        self.model_hash=hashlib.sha256(json.dumps(self.to_dict(include_hash=False),sort_keys=True).encode()).hexdigest(); return self

    def _design(self,z,delta_t,treatment):
        if self.coefficients is None or str(treatment) not in self.coefficients: raise ValueError("unfitted or unseen treatment")
        x=np.asarray(z,float); d=np.asarray(delta_t,float)
        if d.ndim!=0 or not np.isfinite(d) or d<=0: raise ValueError("delta_t must be one positive finite scalar")
        if x.ndim==1: x=x[None,:]
        if x.ndim!=2 or x.shape[1]!=len(self.feature_names): raise ValueError("state dimension mismatch")
        return np.column_stack([np.ones(len(x)),x,np.full(len(x),d)])

    def predict_mean(self,current_state,treatment,delta_t):
        return self._design(current_state,delta_t,treatment) @ np.asarray(self.coefficients[str(treatment)])
    def predict_distribution(self,current_state,treatment,delta_t):
        return {"mean":self.predict_mean(current_state,treatment,delta_t),"covariance":np.asarray(self.residual_covariance[str(treatment)]),"model_hash":self.model_hash,"research_only":True}
    def sample_next(self,current_state,treatment,delta_t,n_samples=1,seed=0):
        if n_samples<1: raise ValueError("n_samples must be positive")
        d=self.predict_distribution(current_state,treatment,delta_t); rng=np.random.default_rng(seed)
        return rng.multivariate_normal(d["mean"][0],d["covariance"],size=n_samples)
    def log_likelihood(self,current_state,next_state,treatment,delta_t):
        from scipy.stats import multivariate_normal
        d=self.predict_distribution(current_state,treatment,delta_t); y=np.asarray(next_state,float)
        return float(np.sum(multivariate_normal.logpdf(y,mean=d["mean"][0],cov=d["covariance"],allow_singular=False)))
    def to_dict(self,include_hash=True):
        d={"version":self.version,"coefficients":self.coefficients,"residual_covariance":self.residual_covariance,"feature_names":self.feature_names,"treatment_levels":self.treatment_levels,"training_metadata":self.training_metadata}
        if include_hash:d["model_hash"]=self.model_hash
        return d
