"""Explicit measurement uncertainty for latent scores."""
from __future__ import annotations
import math
from dataclasses import dataclass
from typing import Tuple

@dataclass(frozen=True)
class UncertaintyComponents:
    measurement_uncertainty: float
    sampling_uncertainty: float
    coverage_penalty: float
    extrapolation_flag: bool
    missing_required_modalities: Tuple[str,...]
    overall_quality: str
    uncertainty_score: float


def uncertainty_score(base_se: float, n_cells: int, gene_coverage: float,
                       modality_coverage: float, reference_cells: int = 500) -> float:
    """Heuristic score; not a calibrated statistical error."""
    if n_cells < 1 or not 0 < gene_coverage <= 1 or not 0 < modality_coverage <= 1:
        return float("inf")
    cell_penalty = math.sqrt(max(1.0, reference_cells / n_cells))
    return float(base_se * cell_penalty / (gene_coverage * modality_coverage))

coverage_adjusted_uncertainty = uncertainty_score  # backwards-compatible alias

def uncertainty_components(base_se, n_cells, gene_coverage, modality_coverage, *, extrapolated=False):
    sampling = float(base_se / max(1.0, n_cells) ** .5)
    penalty = float(1 / max(gene_coverage, 1e-12) / max(modality_coverage, 1e-12))
    return {"measurement_uncertainty": float(base_se), "coverage_penalty": penalty,
            "sampling_uncertainty": sampling, "extrapolation_flag": bool(extrapolated),
            "uncertainty_score": float(base_se * penalty * max(1., (500 / max(n_cells,1)) ** .5))}

def build_uncertainty(base_se,n_cells,gene_coverage,modality_coverage,missing_required_modalities=(),maximum=2.0):
    parts=uncertainty_components(base_se,n_cells,gene_coverage,modality_coverage,extrapolated=bool(missing_required_modalities))
    score=parts["uncertainty_score"]
    if n_cells<1 or gene_coverage<=0 or modality_coverage<=0 or missing_required_modalities: quality="INSUFFICIENT"
    elif score>maximum: quality="LOW_CONFIDENCE"
    else: quality="PASS"
    return UncertaintyComponents(parts["measurement_uncertainty"],parts["sampling_uncertainty"],parts["coverage_penalty"],parts["extrapolation_flag"],tuple(missing_required_modalities),quality,score)


def quality_label(value: float, maximum: float) -> str:
    if not math.isfinite(value): return "INSUFFICIENT"
    return "PASS" if value <= maximum else "LOW_CONFIDENCE"
