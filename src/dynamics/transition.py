"""Treatment-conditioned Gaussian transition model."""
from dataclasses import dataclass
import numpy as np
@dataclass
class GaussianTransition:
    effects: dict[str,np.ndarray] | None=None
    def fit(self,z0,z1,treatment):
        z0,z1=np.asarray(z0,float),np.asarray(z1,float); tr=np.asarray(treatment)
        if len(z0)!=len(z1) or len(z0)!=len(tr): raise ValueError("unaligned transitions")
        self.effects={str(t):(z1[tr==t]-z0[tr==t]).mean(0) for t in np.unique(tr)}; return self
    def predict(self,z,treatment):
        if self.effects is None or treatment not in self.effects: raise ValueError("unfitted or unseen treatment")
        return np.asarray(z,float)+self.effects[treatment]
