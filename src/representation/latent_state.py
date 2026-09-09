"""Typed continuous latent state and tumour-level distribution summaries."""
from __future__ import annotations
from dataclasses import asdict, dataclass, field
from typing import Mapping
import numpy as np


@dataclass(frozen=True)
class LatentFeature:
    value: float
    uncertainty: float
    n_cells: int
    gene_coverage: float
    modality_coverage: float
    qc: str
    measurement_uncertainty: float = float("nan")
    coverage_penalty: float = float("nan")
    sampling_uncertainty: float = float("nan")
    extrapolation_flag: bool = False


@dataclass(frozen=True)
class LatentState:
    sample_id: str
    patient_id: str
    cohort_id: str
    timepoint: str
    features: Mapping[str, LatentFeature]
    available_modalities: tuple[str, ...]
    missing_modalities: tuple[str, ...]
    missing_required_modalities: tuple[str, ...]
    encoder_version: str
    encoder_hash: str

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class TumourStateDistribution:
    """T[t] = p(z | tumour, t), retaining heterogeneity rather than only means."""
    mean: list[float]
    variance: list[float]
    quantiles: Mapping[str, list[float]]
    state_proportions: Mapping[str, float]
    rare_state_abundance: float
    entropy: float
    mixture_components: list[dict] = field(default_factory=list)

    @classmethod
    def from_states(cls, values: np.ndarray, labels: list[str] | None = None,
                    rare_threshold: float = 0.05) -> "TumourStateDistribution":
        x = np.asarray(values, float)
        if x.ndim != 2 or not len(x): raise ValueError("values must be a non-empty cells x dimensions matrix")
        if labels is None: labels = ["unassigned"] * len(x)
        unique, counts = np.unique(labels, return_counts=True)
        props = {str(k): float(v / len(labels)) for k, v in zip(unique, counts)}
        p = np.asarray(list(props.values()))
        return cls(x.mean(0).tolist(), x.var(0, ddof=1).tolist(),
                   {str(q): np.quantile(x, q, axis=0).tolist() for q in (.05, .25, .5, .75, .95)},
                   props, float(p[p < rare_threshold].sum()), float(-(p * np.log(p)).sum()), [])

    @classmethod
    def from_continuous(cls, values: np.ndarray, rare_region_threshold: float = 0.05):
        """Summarise continuous z without treating arbitrary labels as biology.

        Rare regions are defined only as empirical low-density observations and are
        reported as a sampling descriptor, never as a named biological state.
        """
        x = np.asarray(values, float)
        if x.ndim != 2 or len(x) < 2: raise ValueError("at least two continuous states are required")
        from scipy.spatial import distance_matrix
        d = distance_matrix(x, x); np.fill_diagonal(d, np.inf)
        local = np.min(d, axis=1); cutoff = np.quantile(local, 1 - rare_region_threshold)
        p = np.full(len(x), 1 / len(x))
        return cls(x.mean(0).tolist(), x.var(0, ddof=1).tolist(),
                   {str(q): np.quantile(x, q, axis=0).tolist() for q in (.05,.25,.5,.75,.95)},
                   {}, float(np.mean(local >= cutoff)), float(-(p*np.log(p)).sum()),
                   [{"kind":"continuous_low_density_region","threshold":float(cutoff),"n_observations":int(np.sum(local>=cutoff))}])
