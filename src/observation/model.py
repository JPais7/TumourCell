"""Linear-Gaussian p(observation | latent_state, modality)."""
from dataclasses import dataclass
import numpy as np
@dataclass
class LinearGaussianObservation:
    modality: str; loading: np.ndarray; noise_variance: np.ndarray
    def expected(self,latent): return np.asarray(latent) @ np.asarray(self.loading).T
    def reconstruction_error(self,latent,observed): return float(np.mean((self.expected(latent)-observed)**2))
