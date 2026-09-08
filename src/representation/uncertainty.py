"""Explicit measurement uncertainty for latent scores."""
from __future__ import annotations
import math


def coverage_adjusted_uncertainty(base_se: float, n_cells: int, gene_coverage: float,
                                  modality_coverage: float, reference_cells: int = 500) -> float:
    """Inflate uncertainty monotonically for low cell/gene/modality coverage."""
    if n_cells < 1 or not 0 < gene_coverage <= 1 or not 0 < modality_coverage <= 1:
        return float("inf")
    cell_penalty = math.sqrt(max(1.0, reference_cells / n_cells))
    return float(base_se * cell_penalty / (gene_coverage * modality_coverage))


def quality_label(value: float, maximum: float) -> str:
    if not math.isfinite(value): return "INSUFFICIENT"
    return "PASS" if value <= maximum else "LOW_CONFIDENCE"
