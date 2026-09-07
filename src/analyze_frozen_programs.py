#!/usr/bin/env python3
"""Score frozen programs, then analyze response trajectories at patient level."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
FROZEN_NPZ = ROOT / "results/phase1/blind_discovery/program_definitions_v1.npz"
FROZEN_JSON = ROOT / "results/phase1/blind_discovery/program_definitions_v1.json"
PSEUDOBULK = ROOT / "data/processed/GSE246613/malignant_pseudobulk_counts.npz"
OUTDIR = ROOT / "results/phase1/response_analysis"
FIGDIR = ROOT / "figures/phase1"
THRESHOLDS = [25, 50, 100]
PRIMARY_THRESHOLD = 50
TIMES = ["Base", "PD1", "RTPD1"]
STAGES = {"pembrolizumab": ("Base", "PD1"), "radiotherapy_addition": ("PD1", "RTPD1")}
CONTRASTS = {
    "R1_minus_NR": ("R1", "NR"),
    "R2_minus_NR": ("R2", "NR"),
    "R2_minus_R1": ("R2", "R1"),
}
BOOTSTRAPS = 5000
RNG_SEED = 20260906


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def bootstrap_mean_ci(values: np.ndarray, rng: np.random.Generator) -> list[float]:
    if not len(values):
        return [float("nan"), float("nan")]
    draws = rng.choice(values, size=(BOOTSTRAPS, len(values)), replace=True).mean(axis=1)
    return [float(x) for x in np.quantile(draws, [0.025, 0.975])]


def bootstrap_contrast_ci(a: np.ndarray, b: np.ndarray, rng: np.random.Generator) -> list[float]:
    if not len(a) or not len(b):
        return [float("nan"), float("nan")]
    draws_a = rng.choice(a, size=(BOOTSTRAPS, len(a)), replace=True).mean(axis=1)
    draws_b = rng.choice(b, size=(BOOTSTRAPS, len(b)), replace=True).mean(axis=1)
    return [float(x) for x in np.quantile(draws_a - draws_b, [0.025, 0.975])]


def patient_matrices(scores: np.ndarray, patients: np.ndarray, times: np.ndarray, responses: np.ndarray, cells: np.ndarray, threshold: int):
    lookup = {(p, t): i for i, (p, t) in enumerate(zip(patients, times))}
    eligible = sorted(
        p for p in set(patients)
        if all((p, t) in lookup and cells[lookup[(p, t)]] >= threshold for t in TIMES)
    )
    score_cube = np.stack([[scores[lookup[(p, t)]] for t in TIMES] for p in eligible])
    response = np.asarray([responses[lookup[(p, "Base")]] for p in eligible])
    cell_cube = np.asarray([[cells[lookup[(p, t)]] for t in TIMES] for p in eligible])
    return np.asarray(eligible), response, score_cube, cell_cube


def contrast_record(values: np.ndarray, response: np.ndarray, group_a: str, group_b: str, rng: np.random.Generator):
    a, b = values[response == group_a], values[response == group_b]
    effect = float(a.mean() - b.mean())
    if len(a) >= 2 and len(b) >= 2:
        pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
        standardized = float(effect / pooled) if pooled > 0 else None
    else:
        standardized = None
    return {
        "group_a": group_a,
        "group_b": group_b,
        "n_a": int(len(a)),
        "n_b": int(len(b)),
        "mean_difference": effect,
        "bootstrap_95_ci": bootstrap_contrast_ci(a, b, rng),
        "standardized_mean_difference": standardized,
    }


def leave_one_out(values: np.ndarray, response: np.ndarray, patients: np.ndarray, group_a: str, group_b: str):
    full = values[response == group_a].mean() - values[response == group_b].mean()
    records = []
    for patient in patients[(response == group_a) | (response == group_b)]:
        keep = patients != patient
        a, b = values[keep & (response == group_a)], values[keep & (response == group_b)]
        if not len(a) or not len(b):
            continue
        effect = float(a.mean() - b.mean())
        records.append((str(patient), effect, abs(effect - full)))
    same_sign = np.mean([np.sign(effect) == np.sign(full) for _, effect, _ in records]) if records else np.nan
    influential = sorted(records, key=lambda item: item[2], reverse=True)[:3]
    return {
        "full_effect": float(full),
        "leave_one_out_direction_stability": float(same_sign),
        "most_influential": [
            {"patient": patient, "leave_one_out_effect": effect, "absolute_effect_change": change}
            for patient, effect, change in influential
        ],
    }


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    manifest = json.loads(FROZEN_JSON.read_text())
    if not manifest["frozen"] or sha256(FROZEN_NPZ) != manifest["frozen_npz_sha256"]:
        raise RuntimeError("frozen program definition failed integrity check")

    frozen = np.load(FROZEN_NPZ)
    pb = np.load(PSEUDOBULK)
    genes = frozen["genes"].astype(str)
    all_genes = pb["genes"].astype(str)
    lookup = {gene: i for i, gene in enumerate(all_genes)}
    columns = np.asarray([lookup[gene] for gene in genes])
    counts = pb["counts"][:, columns].astype(np.float64)
    library = pb["counts"].sum(axis=1, keepdims=True)
    log_cpm = np.log1p(counts / np.maximum(library, 1) * 1_000_000)
    scaled = log_cpm / frozen["gene_scale"]
    raw_scores = scaled @ frozen["weights"].T
    score_mean = raw_scores.mean(axis=0)
    score_sd = raw_scores.std(axis=0, ddof=1)
    scores = (raw_scores - score_mean) / np.maximum(score_sd, 1e-12)
    program_ids = np.asarray([f"P{i}" for i in range(1, scores.shape[1] + 1)])
    labels = frozen["labels"].astype(str)
    patients = pb["patient"].astype(str)
    times = pb["treatment"].astype(str)
    responses = pb["response_group"].astype(str)  # First response-variable access occurs here.
    cells = pb["malignant_cells"].astype(int)

    np.savez_compressed(
        OUTDIR / "program_scores_v1.npz",
        patient=patients,
        treatment=times,
        response_group=responses,
        malignant_cells=cells,
        program_ids=program_ids,
        labels=labels,
        scores=scores,
        raw_scores=raw_scores,
        score_center=score_mean,
        score_scale=score_sd,
        frozen_definition_sha256=np.asarray(manifest["frozen_npz_sha256"]),
    )

    rng = np.random.default_rng(RNG_SEED)
    analysis = {
        "frozen_definition_sha256": manifest["frozen_npz_sha256"],
        "response_opened_after_freeze": True,
        "primary_threshold": PRIMARY_THRESHOLD,
        "bootstrap_replicates": BOOTSTRAPS,
        "bootstrap_seed": RNG_SEED,
        "thresholds": {},
    }
    primary_payload = None
    for threshold in THRESHOLDS:
        pt, response, cube, cell_cube = patient_matrices(scores, patients, times, responses, cells, threshold)
        deltas = {
            stage: cube[:, TIMES.index(end), :] - cube[:, TIMES.index(start), :]
            for stage, (start, end) in STAGES.items()
        }
        payload = {
            "patients": pt.tolist(),
            "response_counts": dict(Counter(response)),
            "programs": {},
        }
        for j, (program, label) in enumerate(zip(program_ids, labels)):
            program_result = {"label": label, "stages": {}}
            for stage, matrix in deltas.items():
                values = matrix[:, j]
                groups = {}
                for group in ["R1", "R2", "NR"]:
                    selected = values[response == group]
                    groups[group] = {
                        "n": int(len(selected)),
                        "mean_delta": float(selected.mean()),
                        "median_delta": float(np.median(selected)),
                        "bootstrap_mean_95_ci": bootstrap_mean_ci(selected, rng),
                    }
                contrasts = {
                    name: contrast_record(values, response, a, b, rng)
                    for name, (a, b) in CONTRASTS.items()
                }
                influence = {
                    name: leave_one_out(values, response, pt, a, b)
                    for name, (a, b) in CONTRASTS.items()
                }
                depth = np.minimum(
                    cell_cube[:, TIMES.index(STAGES[stage][0])],
                    cell_cube[:, TIMES.index(STAGES[stage][1])],
                )
                rho, pvalue = spearmanr(values, depth)
                program_result["stages"][stage] = {
                    "groups": groups,
                    "contrasts": contrasts,
                    "leave_one_patient_out": influence,
                    "depth_spearman_rho": float(rho),
                    "depth_spearman_p": float(pvalue),
                }
            payload["programs"][str(program)] = program_result
        analysis["thresholds"][str(threshold)] = payload
        if threshold == PRIMARY_THRESHOLD:
            primary_payload = (pt, response, cube, deltas)

    # Residualize non-confounder deltas against frozen stress (P2) and cycle (P3) deltas.
    pt, response, cube, deltas = primary_payload
    adjusted = {}
    for stage, matrix in deltas.items():
        design = np.column_stack([np.ones(len(pt)), matrix[:, 1], matrix[:, 2]])
        adjusted[stage] = {}
        for j, program in enumerate(program_ids):
            if j in {1, 2}:
                continue
            residual = matrix[:, j] - design @ np.linalg.lstsq(design, matrix[:, j], rcond=None)[0]
            adjusted[stage][str(program)] = {
                name: contrast_record(residual, response, a, b, rng)
                for name, (a, b) in CONTRASTS.items()
            }
    analysis["primary_adjusted_for_P2_stress_and_P3_cycle"] = adjusted
    (OUTDIR / "patient_level_contrasts_v1.json").write_text(json.dumps(analysis, indent=2) + "\n")

    # Full heatmap: no result-driven program omission.
    fig, axes = plt.subplots(1, 2, figsize=(10, 5), sharey=True)
    for ax, (stage, matrix) in zip(axes, deltas.items()):
        means = np.asarray([[matrix[response == group, j].mean() for group in ["R1", "R2", "NR"]] for j in range(len(program_ids))])
        image = ax.imshow(means, cmap="coolwarm", vmin=-1.5, vmax=1.5, aspect="auto")
        ax.set_xticks(range(3), ["R1", "R2", "NR"])
        ax.set_yticks(range(len(program_ids)), [f"{p} {label}" for p, label in zip(program_ids, labels)])
        ax.set_title(stage.replace("_", " "))
        for i in range(means.shape[0]):
            for j in range(means.shape[1]):
                ax.text(j, i, f"{means[i,j]:.2f}", ha="center", va="center", fontsize=7)
    fig.subplots_adjust(left=0.30, right=0.84, wspace=0.15)
    color_axis = fig.add_axes([0.87, 0.22, 0.025, 0.60])
    fig.colorbar(image, cax=color_axis, label="Mean within-patient change (SD units)")
    fig.savefig(FIGDIR / "program_longitudinal_changes_primary.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIGDIR / "program_longitudinal_changes_primary.pdf", bbox_inches="tight")
    print(json.dumps({"primary_patients": len(pt), "response_counts": dict(Counter(response))}, indent=2))


if __name__ == "__main__":
    main()
