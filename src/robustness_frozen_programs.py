#!/usr/bin/env python3
"""Sensitivity analysis for frozen programs without redefining the primary model."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "results/phase1/blind_discovery/program_definitions_v1.npz"
CANDIDATES = ROOT / "results/phase1/blind_discovery/candidate_program_weights.npz"
PSEUDOBULK = ROOT / "data/processed/GSE246613/malignant_pseudobulk_counts.npz"
PRIMARY = ROOT / "results/phase1/response_analysis/patient_level_contrasts_v1.json"
GENE_SENSITIVITY = ROOT / "results/phase1/blind_discovery/gene_count_sensitivity_weights.npz"
OUTPUT = ROOT / "results/phase1/response_analysis/robustness_v1.json"
TIMES = ["Base", "PD1", "RTPD1"]


def normalized_rows(matrix: np.ndarray) -> np.ndarray:
    return matrix / np.maximum(matrix.sum(axis=1, keepdims=True), 1e-12)


def cosine_rows(matrix: np.ndarray) -> np.ndarray:
    return matrix / np.maximum(np.linalg.norm(matrix, axis=1, keepdims=True), 1e-12)


def aligned_weights(reference: np.ndarray, alternative: np.ndarray):
    similarity = cosine_rows(reference) @ cosine_rows(alternative).T
    rows, cols = linear_sum_assignment(-similarity)
    mapping = {int(row): int(col) for row, col in zip(rows, cols)}
    return mapping, {f"P{row + 1}": float(similarity[row, col]) for row, col in zip(rows, cols)}


def score(counts: np.ndarray, library: np.ndarray, scale: np.ndarray, weights: np.ndarray, normalization: str):
    if normalization == "log_cpm":
        expression = np.log1p(counts / np.maximum(library, 1) * 1_000_000)
    elif normalization == "log_tp10k":
        expression = np.log1p(counts / np.maximum(library, 1) * 10_000)
    elif normalization == "sqrt_cpm":
        expression = np.sqrt(counts / np.maximum(library, 1) * 1_000_000)
    else:
        raise ValueError(normalization)
    values = expression / scale @ normalized_rows(weights).T
    return (values - values.mean(axis=0)) / np.maximum(values.std(axis=0, ddof=1), 1e-12)


def effect(scores, patients, times, responses, cells, threshold, stage, program):
    lookup = {(p, t): i for i, (p, t) in enumerate(zip(patients, times))}
    eligible = sorted(p for p in set(patients) if all((p, t) in lookup and cells[lookup[(p, t)]] >= threshold for t in TIMES))
    start, end = {"pembrolizumab": ("Base", "PD1"), "radiotherapy_addition": ("PD1", "RTPD1")}[stage]
    delta = np.asarray([scores[lookup[(p, end)], program] - scores[lookup[(p, start)], program] for p in eligible])
    group = np.asarray([responses[lookup[(p, "Base")]] for p in eligible])
    a, b = delta[group == "R2"], delta[group == "NR"]
    return float(a.mean() - b.mean()), len(a), len(b)


def main() -> None:
    frozen = np.load(FROZEN)
    candidates = np.load(CANDIDATES)
    pb = np.load(PSEUDOBULK)
    genes = frozen["genes"].astype(str)
    all_genes = pb["genes"].astype(str)
    columns = np.asarray([{gene: i for i, gene in enumerate(all_genes)}[gene] for gene in genes])
    counts = pb["counts"][:, columns]
    library = pb["counts"].sum(axis=1, keepdims=True)
    patients = pb["patient"].astype(str)
    times = pb["treatment"].astype(str)
    responses = pb["response_group"].astype(str)
    cells = pb["malignant_cells"].astype(int)
    reference = frozen["weights"]

    models = {"rank9_seed11_primary": reference}
    for seed in [29, 47, 71, 101]:
        models[f"rank9_seed{seed}"] = candidates[f"rank_9_seed_{seed}"]
    for rank in [8, 10]:
        for seed in [11, 29, 47, 71, 101]:
            models[f"rank{rank}_seed{seed}"] = candidates[f"rank_{rank}_seed_{seed}"]

    records = []
    for name, weights in models.items():
        mapping, similarities = aligned_weights(reference, weights)
        model_scores = score(counts, library, frozen["gene_scale"], weights, "log_cpm")
        for ref_program, alt_program in mapping.items():
            for threshold in [25, 50, 100]:
                for stage in ["pembrolizumab", "radiotherapy_addition"]:
                    value, n_r2, n_nr = effect(model_scores, patients, times, responses, cells, threshold, stage, alt_program)
                    records.append({
                        "variant": name,
                        "reference_program": f"P{ref_program + 1}",
                        "matched_alternative_program": int(alt_program + 1),
                        "component_cosine": similarities[f"P{ref_program + 1}"],
                        "threshold": threshold,
                        "stage": stage,
                        "contrast": "R2_minus_NR",
                        "effect": value,
                        "n_R2": n_r2,
                        "n_NR": n_nr,
                    })

    normalization_records = []
    for normalization in ["log_cpm", "log_tp10k", "sqrt_cpm"]:
        model_scores = score(counts, library, frozen["gene_scale"], reference, normalization)
        for program in range(len(reference)):
            for threshold in [25, 50, 100]:
                for stage in ["pembrolizumab", "radiotherapy_addition"]:
                    value, n_r2, n_nr = effect(model_scores, patients, times, responses, cells, threshold, stage, program)
                    normalization_records.append({
                        "normalization": normalization,
                        "program": f"P{program + 1}",
                        "threshold": threshold,
                        "stage": stage,
                        "contrast": "R2_minus_NR",
                        "effect": value,
                        "n_R2": n_r2,
                        "n_NR": n_nr,
                    })

    gene_count_records = []
    sensitivity = np.load(GENE_SENSITIVITY)
    frozen_gene_lookup = {gene: i for i, gene in enumerate(genes)}
    all_gene_lookup = {gene: i for i, gene in enumerate(all_genes)}
    for gene_count in [1500, 2500]:
        alt_genes = sensitivity[f"genes_{gene_count}"].astype(str)
        alt_weights = sensitivity[f"weights_{gene_count}"]
        alt_scale = sensitivity[f"gene_scale_{gene_count}"]
        common = sorted(set(genes) & set(alt_genes))
        ref_common = reference[:, [frozen_gene_lookup[g] for g in common]]
        alt_lookup = {gene: i for i, gene in enumerate(alt_genes)}
        alt_common = alt_weights[:, [alt_lookup[g] for g in common]]
        mapping, similarities = aligned_weights(ref_common, alt_common)
        alt_columns = np.asarray([all_gene_lookup[g] for g in alt_genes])
        alt_scores = score(pb["counts"][:, alt_columns], library, alt_scale, alt_weights, "log_cpm")
        for ref_program, alt_program in mapping.items():
            for threshold in [25, 50, 100]:
                for stage in ["pembrolizumab", "radiotherapy_addition"]:
                    value, n_r2, n_nr = effect(alt_scores, patients, times, responses, cells, threshold, stage, alt_program)
                    gene_count_records.append({
                        "selected_genes": gene_count,
                        "common_genes_with_primary": len(common),
                        "reference_program": f"P{ref_program + 1}",
                        "matched_alternative_program": int(alt_program + 1),
                        "component_cosine_on_common_genes": similarities[f"P{ref_program + 1}"],
                        "threshold": threshold,
                        "stage": stage,
                        "contrast": "R2_minus_NR",
                        "effect": value,
                        "n_R2": n_r2,
                        "n_NR": n_nr,
                    })

    primary = json.loads(PRIMARY.read_text())
    adjusted = primary["primary_adjusted_for_P2_stress_and_P3_cycle"]
    summary = {}
    for program in [f"P{i}" for i in range(1, 10)]:
        relevant = [x for x in records if x["reference_program"] == program and x["threshold"] == 50 and x["stage"] == "radiotherapy_addition"]
        norm = [x for x in normalization_records if x["program"] == program and x["threshold"] == 50 and x["stage"] == "radiotherapy_addition"]
        primary_effect = next(x["effect"] for x in relevant if x["variant"] == "rank9_seed11_primary")
        summary[program] = {
            "primary_R2_minus_NR_radiotherapy_effect": primary_effect,
            "rank_seed_variants_same_direction_fraction": float(np.mean([np.sign(x["effect"]) == np.sign(primary_effect) for x in relevant])),
            "normalizations_same_direction_fraction": float(np.mean([np.sign(x["effect"]) == np.sign(primary_effect) for x in norm])),
            "gene_counts_same_direction_fraction": float(np.mean([
                np.sign(x["effect"]) == np.sign(primary_effect)
                for x in gene_count_records
                if x["reference_program"] == program and x["threshold"] == 50 and x["stage"] == "radiotherapy_addition"
            ])),
            "adjusted_effect": adjusted["radiotherapy_addition"].get(program, {}).get("R2_minus_NR", {}).get("mean_difference"),
        }
    payload = {
        "primary_model_unchanged": True,
        "contrast_scope": "R2 versus NR; R1 excluded from inferential sensitivity because primary n=1",
        "rank_seed_alignment": "Hungarian matching of component cosine similarity to frozen rank-9 seed-11 programs",
        "model_records": records,
        "normalization_records": normalization_records,
        "gene_count_records": gene_count_records,
        "program_summary": summary,
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
