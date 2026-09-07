#!/usr/bin/env python3
"""Apply the frozen Phase-1 P8 definition to ARTEMIS without retuning."""

from __future__ import annotations

import hashlib
import json
import platform
from collections import Counter
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.sparse import csr_matrix


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/ARTEMIS/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad"
SPEC_FILE = ROOT / "results/phase2/p8_definition_frozen.json"
PROTOCOL_FILE = ROOT / "docs/phase2_validation_protocol.md"
OUT_JSON = ROOT / "results/phase2/artemis_p8_validation.json"
OUT_NPZ = ROOT / "results/phase2/artemis_p8_patient_scores.npz"
FIGURE = ROOT / "results/figures/phase2_p8_artemis.png"
EXPECTED_INPUT_SHA256 = "d40a2b210ddf8fe0bd46118900412493c6d11d3429ffae4c311f7469f7d10da0"
EXPECTED_SPEC_SHA256 = "185cde59ee15f1fac893b8367cfd1320dd6e1559d8e4894ab060cebfc4ed744e"
EXPECTED_PROTOCOL_SHA256 = "583e704151fc3f17536e9b9b9912f3922391aab9ef4d6a3afbbded3087350897"
SEED = 20260906


def sha256(path: Path, chunk: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while data := handle.read(chunk):
            digest.update(data)
    return digest.hexdigest()


def decode(values) -> np.ndarray:
    return np.asarray([
        value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def obs_column(obs: h5py.Group, name: str) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = decode(obj["categories"][:])
        codes = obj["codes"][:]
        return np.asarray([categories[code] if code >= 0 else "<NA>" for code in codes])
    return decode(obj[:])


def aggregate_tumour_counts() -> dict:
    with h5py.File(INPUT, "r") as handle:
        obs = handle["obs"]
        author_type = obs_column(obs, "author_cell_type")
        patient_cell = obs_column(obs, "donor_id")
        outcome_cell = obs_column(obs, "pCR_status")
        assay_cell = obs_column(obs, "assay")
        tumour = author_type == "Tumor"
        patients = np.unique(patient_cell[tumour])
        patient_index = {name: i for i, name in enumerate(patients)}
        cell_group = np.full(len(tumour), -1, dtype=np.int16)
        for name, index in patient_index.items():
            cell_group[tumour & (patient_cell == name)] = index

        genes = decode(handle["raw/var/gene_symbols"][:])
        counts = np.zeros((len(patients), len(genes)), dtype=np.float64)
        matrix = handle["raw/X"]
        indptr = matrix["indptr"][:]
        data_ds, indices_ds = matrix["data"], matrix["indices"]
        for start in range(0, len(tumour), 2048):
            end = min(start + 2048, len(tumour))
            chosen = cell_group[start:end] >= 0
            if not chosen.any():
                continue
            left, right = int(indptr[start]), int(indptr[end])
            chunk = csr_matrix(
                (data_ds[left:right], indices_ds[left:right], indptr[start : end + 1] - left),
                shape=(end - start, len(genes)),
            )
            local = cell_group[start:end]
            for index in np.unique(local[chosen]):
                counts[index] += np.asarray(chunk[local == index].sum(axis=0)).ravel()

        outcomes, assays, n_cells = [], [], []
        for name in patients:
            mask = tumour & (patient_cell == name)
            found_outcome = sorted(set(outcome_cell[mask]))
            found_assay = sorted(set(assay_cell[mask]))
            if len(found_outcome) != 1 or len(found_assay) != 1:
                raise ValueError(f"Inconsistent metadata for {name}: {found_outcome}, {found_assay}")
            outcomes.append(found_outcome[0])
            assays.append(found_assay[0])
            n_cells.append(int(mask.sum()))

    return {
        "counts": counts,
        "genes": genes,
        "patient": patients,
        "outcome": np.asarray(outcomes),
        "assay": np.asarray(assays),
        "n_cells": np.asarray(n_cells),
    }


def p8_score(counts: np.ndarray, library: np.ndarray, spec: dict, gene_index: dict, mode: str) -> np.ndarray:
    source_genes = np.asarray(spec["genes"])
    present = np.asarray([gene in gene_index for gene in source_genes])
    source_cols = np.flatnonzero(present)
    target_cols = np.asarray([gene_index[gene] for gene in source_genes[present]])
    values = counts[:, target_cols]
    if mode == "log_cpm":
        transformed = np.log1p(values / library[:, None] * 1_000_000.0)
    elif mode == "log_tp10k":
        transformed = np.log1p(values / library[:, None] * 10_000.0)
    elif mode == "sqrt_cpm":
        transformed = np.sqrt(values / library[:, None] * 1_000_000.0)
    else:
        raise ValueError(mode)
    weights = np.asarray(spec["weights"])[source_cols]
    weights = weights / weights.sum()
    scales = np.asarray(spec["gene_scale"])[source_cols]
    return (transformed / scales) @ weights


def zscore_on_technical(raw_score: np.ndarray, technical: np.ndarray) -> np.ndarray:
    mean = raw_score[technical].mean()
    sd = raw_score[technical].std(ddof=1)
    return (raw_score - mean) / sd


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    n1, n0 = len(a), len(b)
    pooled = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n0 - 1) * b.var(ddof=1)) / (n1 + n0 - 2))
    d = (a.mean() - b.mean()) / pooled
    return float(d * (1 - 3 / (4 * (n1 + n0) - 9)))


def primary_stats(score: np.ndarray, eligible: np.ndarray, outcome: np.ndarray, patient: np.ndarray) -> dict:
    use = eligible & np.isin(outcome, ["pCR", "RD"])
    a = score[use & (outcome == "pCR")]
    b = score[use & (outcome == "RD")]
    effect = float(a.mean() - b.mean())
    rng = np.random.default_rng(SEED)
    boot = np.empty(10_000)
    for i in range(len(boot)):
        boot[i] = rng.choice(a, len(a), replace=True).mean() - rng.choice(b, len(b), replace=True).mean()
    loo = []
    for idx in np.flatnonzero(use):
        keep = use.copy()
        keep[idx] = False
        value = float(score[keep & (outcome == "pCR")].mean() - score[keep & (outcome == "RD")].mean())
        loo.append({"patient": str(patient[idx]), "outcome": str(outcome[idx]), "effect": value,
                    "absolute_change": abs(value - effect)})
    loo.sort(key=lambda item: item["absolute_change"], reverse=True)
    return {
        "n": int(use.sum()), "n_pCR": int(len(a)), "n_RD": int(len(b)),
        "mean_pCR": float(a.mean()), "mean_RD": float(b.mean()),
        "mean_difference_pCR_minus_RD": effect,
        "bootstrap_ci95": [float(x) for x in np.quantile(boot, [0.025, 0.975])],
        "median_difference_pCR_minus_RD": float(np.median(a) - np.median(b)),
        "hedges_g": hedges_g(a, b),
        "sd_pCR": float(a.std(ddof=1)), "sd_RD": float(b.std(ddof=1)),
        "loo_positive_fraction": float(np.mean([item["effect"] > 0 for item in loo])),
        "loo_effect_range": [float(min(item["effect"] for item in loo)), float(max(item["effect"] for item in loo))],
        "most_influential": loo[:3],
        "eligible_patient_ids": [str(x) for x in patient[use]],
    }


def adjusted_effect(score, eligible, outcome, library, n_cells, assay) -> float:
    use = eligible & np.isin(outcome, ["pCR", "RD"])
    y = score[use]
    pcr = (outcome[use] == "pCR").astype(float)
    cov = np.column_stack([
        np.ones(use.sum()), pcr,
        (np.log(library[use]) - np.log(library[use]).mean()) / np.log(library[use]).std(ddof=1),
        (n_cells[use] - n_cells[use].mean()) / n_cells[use].std(ddof=1),
        (assay[use] == "10x 3' v2").astype(float),
    ])
    return float(np.linalg.lstsq(cov, y, rcond=None)[0][1])


def main() -> None:
    hashes = {"input": sha256(INPUT), "p8_spec": sha256(SPEC_FILE), "protocol": sha256(PROTOCOL_FILE)}
    expected = {"input": EXPECTED_INPUT_SHA256, "p8_spec": EXPECTED_SPEC_SHA256, "protocol": EXPECTED_PROTOCOL_SHA256}
    if hashes != expected:
        raise RuntimeError(f"Frozen/input hash mismatch: {hashes}")
    spec = json.loads(SPEC_FILE.read_text())
    data = aggregate_tumour_counts()
    counts, genes = data["counts"], data["genes"]
    patient, outcome, assay, n_cells = data["patient"], data["outcome"], data["assay"], data["n_cells"]
    library = counts.sum(axis=1)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    present = np.asarray([gene in gene_index for gene in spec["genes"]])
    weights = np.asarray(spec["weights"])
    top50 = np.argsort(-weights)[:50]
    coverage = {
        "genes_present": int(present.sum()), "genes_required": 1600,
        "weight_present": float(weights[present].sum()), "weight_required": 0.80,
        "top50_present": int(present[top50].sum()), "top50_required": 40,
        "missing_top50": [spec["genes"][i] for i in top50 if not present[i]],
    }
    coverage["pass"] = bool(coverage["genes_present"] >= 1600 and coverage["weight_present"] >= .8 and coverage["top50_present"] >= 40)
    if not coverage["pass"]:
        raise RuntimeError(f"Frozen coverage rules failed: {coverage}")

    raw_scores = {mode: p8_score(counts, library, spec, gene_index, mode) for mode in ["log_cpm", "log_tp10k", "sqrt_cpm"]}
    analyses = {}
    standardized = {}
    for threshold in [25, 50, 100]:
        technical = n_cells >= threshold
        score = zscore_on_technical(raw_scores["log_cpm"], technical)
        standardized[threshold] = score
        analyses[f"threshold_{threshold}_log_cpm"] = primary_stats(score, technical, outcome, patient)
    technical50 = n_cells >= 50
    for mode in ["log_tp10k", "sqrt_cpm"]:
        score = zscore_on_technical(raw_scores[mode], technical50)
        analyses[f"threshold_50_{mode}"] = primary_stats(score, technical50, outcome, patient)
    analyses["threshold_50_log_cpm_adjusted_effect"] = adjusted_effect(
        standardized[50], technical50, outcome, library, n_cells, assay
    )
    primary = analyses["threshold_50_log_cpm"]
    sensitivity_effects = [
        analyses["threshold_25_log_cpm"]["mean_difference_pCR_minus_RD"],
        analyses["threshold_100_log_cpm"]["mean_difference_pCR_minus_RD"],
        analyses["threshold_50_log_tp10k"]["mean_difference_pCR_minus_RD"],
        analyses["threshold_50_sqrt_cpm"]["mean_difference_pCR_minus_RD"],
        analyses["threshold_50_log_cpm_adjusted_effect"],
    ]
    if primary["bootstrap_ci95"][0] > 0 and primary["loo_positive_fraction"] >= .90 and all(x > 0 for x in sensitivity_effects):
        verdict = "REPLICAÇÃO PARCIAL"
    elif primary["bootstrap_ci95"][1] < 0 and all(x < 0 for x in sensitivity_effects):
        verdict = "NÃO REPLICAÇÃO"
    else:
        verdict = "INCONCLUSIVO"

    result = {
        "verdict": verdict, "endpoint_comparability": "partial: baseline pCR versus RD, not longitudinal R2 versus NR",
        "hashes": hashes, "dataset_id": "e94bd3cc-6271-424a-baac-12f8eb320a0e",
        "dataset": {
            "cells": 427823, "tumour_cells": int(n_cells.sum()), "patients_total": 101,
            "patients_with_tumour": int(len(patient)), "patient_outcomes_with_tumour": dict(Counter(outcome)),
            "tumour_cells_by_assay": dict(Counter(np.repeat(assay, n_cells))),
            "tumour_cell_threshold_counts": {str(t): int((n_cells >= t).sum()) for t in [25, 50, 100]},
        },
        "coverage": coverage, "primary": primary, "analyses": analyses,
        "sensitivity_all_preserve_primary_direction": bool(all(np.sign(x) == np.sign(primary["mean_difference_pCR_minus_RD"]) for x in sensitivity_effects)),
        "software": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "h5py": h5py.__version__},
        "seed": SEED,
    }
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    np.savez_compressed(OUT_NPZ, patient=patient, outcome=outcome, assay=assay, malignant_cells=n_cells,
                        library_size=library, p8_log_cpm=raw_scores["log_cpm"], p8_z_primary=standardized[50])

    use = technical50 & np.isin(outcome, ["pCR", "RD"])
    fig, ax = plt.subplots(figsize=(5.6, 4.5))
    rng = np.random.default_rng(SEED)
    for x, label, color in [(0, "pCR", "#2a9d8f"), (1, "RD", "#8e5ba6")]:
        vals = standardized[50][use & (outcome == label)]
        ax.boxplot(vals, positions=[x], widths=.5, patch_artist=True, showfliers=False,
                   boxprops={"facecolor": color, "alpha": .35}, medianprops={"color": "black"})
        ax.scatter(x + rng.uniform(-.13, .13, len(vals)), vals, s=18, alpha=.7, color=color, edgecolor="none")
    ax.set_xticks([0, 1], [f"pCR (n={primary['n_pCR']})", f"RD (n={primary['n_RD']})"])
    ax.set_ylabel("P8 congelado (SD ARTEMIS)")
    lo, hi = primary["bootstrap_ci95"]
    ax.set_title(f"pCR−RD = {primary['mean_difference_pCR_minus_RD']:.2f} SD; IC95% {lo:.2f} a {hi:.2f}")
    ax.axhline(0, color="0.75", lw=1, zorder=0)
    fig.tight_layout()
    FIGURE.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)
    print(json.dumps({"verdict": verdict, "primary": primary, "coverage": coverage}, indent=2))


if __name__ == "__main__":
    main()
