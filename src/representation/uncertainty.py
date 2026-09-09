"""Explicit measurement uncertainty for latent scores."""
from __future__ import annotations
import math


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


def quality_label(value: float, maximum: float) -> str:
    if not math.isfinite(value): return "INSUFFICIENT"
    return "PASS" if value <= maximum else "LOW_CONFIDENCE"
