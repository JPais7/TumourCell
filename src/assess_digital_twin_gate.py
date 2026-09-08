#!/usr/bin/env python3
"""Audit whether current evidence permits a patient-specific TNBC simulator."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


ROOT = Path(__file__).resolve().parents[1]
SEED = 20260908


def auc(y: np.ndarray, score: np.ndarray) -> float:
    """AUROC by the Mann-Whitney rank identity; y=1 is favourable outcome."""
    pos, neg = y == 1, y == 0
    ranks = rankdata(score)
    return float((ranks[pos].sum() - pos.sum() * (pos.sum() + 1) / 2) / (pos.sum() * neg.sum()))


def stratified_bootstrap_auc(y: np.ndarray, score: np.ndarray, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    pos, neg = np.flatnonzero(y == 1), np.flatnonzero(y == 0)
    draws = np.asarray([
        auc(
            np.r_[np.ones(len(pos), dtype=int), np.zeros(len(neg), dtype=int)],
            np.r_[score[rng.choice(pos, len(pos), replace=True)], score[rng.choice(neg, len(neg), replace=True)]],
        )
        for _ in range(10000)
    ])
    observed = auc(y, score)
    return {
        "auc": observed,
        "bootstrap_ci95": [float(x) for x in np.quantile(draws, [.025, .975])],
        "n_favourable": int((y == 1).sum()),
        "n_residual": int((y == 0).sum()),
    }


def main() -> None:
    score_path = ROOT / "results/phase3d/gse260693_frozen_scores.csv"
    with score_path.open() as handle:
        rows = [row for row in csv.DictReader(handle) if row["timepoint"] == "pre"]
    y = np.asarray([row["residual_tumor"].lower() == "no" for row in rows], dtype=int)
    discrimination = {}
    for index, signature in enumerate(["State_03", "P8"]):
        values = np.asarray([float(row[signature]) for row in rows])
        result = stratified_bootstrap_auc(y, values, SEED + index)
        result["locked_direction"] = "higher score predicts favourable pathological outcome"
        result["level_2_criterion"] = "external AUROC 95% CI lower bound > 0.50"
        result["passes"] = result["bootstrap_ci95"][0] > 0.5
        discrimination[signature] = result

    gates = {
        "G1_locked_representation_transport": {
            "status": "partial",
            "evidence": "100/100 genes for each state and 97.47% of P8 weight in GSE260693; bulk scores are not malignant-cell-specific",
        },
        "G2_independent_patient_level_prediction": {
            "status": "pass" if any(x["passes"] for x in discrimination.values()) else "fail",
            "evidence": discrimination,
        },
        "G3_treatment_conditioned_transition_model": {
            "status": "blocked",
            "evidence": "NeoTRIP treatment arm and pCR are unavailable; GSE260693 regimens are heterogeneous and has 19 pairs",
        },
        "G4_clone_resolved_dynamics": {
            "status": "blocked",
            "evidence": "no matched clone-resolved longitudinal validation in the analysed cohorts",
        },
        "G5_spatial_external_validation": {
            "status": "blocked",
            "evidence": "no adequate multi-patient whole-transcriptome spatial TNBC response cohort is currently available",
        },
    }
    allowed = gates["G2_independent_patient_level_prediction"]["status"] == "pass"
    payload = {
        "analysis": "digital twin readiness gate v0.1",
        "date": "2026-09-08",
        "unit_of_inference": "patient",
        "outcome": "no residual tumour versus residual tumour in external GSE260693 baseline samples",
        "no_refitting": True,
        "threshold_note": "conservative engineering gate; not retrospectively claimed as preregistered",
        "gates": gates,
        "virtual_perturbation_allowed": allowed,
        "patient_specific_treatment_recommendation_allowed": False,
        "permitted_next_artifact": "observational molecular state card with uncertainty and explicit missingness",
    }
    out = ROOT / "results/phase4"
    out.mkdir(parents=True, exist_ok=True)
    (out / "digital_twin_readiness_gate.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
