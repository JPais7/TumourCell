#!/usr/bin/env python3
"""Response-blind malignant gene selection and NMF rank/seed comparison."""

from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix, save_npz, vstack
from sklearn.decomposition import NMF


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE246613/GSE246613_PembroRT_non_immune_cells.h5ad"
OUTDIR = ROOT / "results/phase1/blind_discovery"
FIGDIR = ROOT / "figures/phase1"
CACHE = ROOT / "data/processed/GSE246613/blind_nmf_input.npz"

N_HVG = 2000
MIN_DETECTION = 0.01
MAX_DETECTION = 0.90
MEAN_BINS = 20
MAX_CELLS_PER_PATIENT_TIME = 100
SAMPLING_SEED = 246613
RANKS = list(range(4, 11))
SEEDS = [11, 29, 47, 71, 101]
CHUNK_ROWS = 2048
MAX_ITER = 500


def decode(values) -> np.ndarray:
    return np.asarray([
        x.decode("utf-8", errors="replace") if isinstance(x, bytes) else str(x)
        for x in values
    ])


def obs_column(obs: h5py.Group, name: str) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = decode(obj["categories"][:])
        codes = obj["codes"][:]
        return np.asarray([categories[code] if code >= 0 else "" for code in codes])
    return decode(obj[:])


def collapse_patient(value: str) -> str:
    return "Patient03" if value in {"Patient03T1", "Patient03T2"} else value


def read_csr_chunk(group: h5py.Group, indptr: np.ndarray, start: int, end: int, n_genes: int) -> csr_matrix:
    left, right = int(indptr[start]), int(indptr[end])
    return csr_matrix(
        (
            group["data"][left:right],
            group["indices"][left:right],
            indptr[start : end + 1] - left,
        ),
        shape=(end - start, n_genes),
    )


def standardized_variance(mean: np.ndarray, variance: np.ndarray, eligible: np.ndarray) -> np.ndarray:
    score = np.full(len(mean), -np.inf)
    idx = np.flatnonzero(eligible)
    order = idx[np.argsort(mean[idx])]
    for members in np.array_split(order, MEAN_BINS):
        values = np.log1p(variance[members])
        spread = values.std(ddof=1)
        score[members] = (values - values.mean()) / spread if spread > 0 else 0.0
    return score


def matched_stability(reference: np.ndarray, other: np.ndarray) -> float:
    ref = reference / np.maximum(np.linalg.norm(reference, axis=1, keepdims=True), 1e-12)
    alt = other / np.maximum(np.linalg.norm(other, axis=1, keepdims=True), 1e-12)
    similarity = ref @ alt.T
    rows, cols = linear_sum_assignment(-similarity)
    return float(np.mean(similarity[rows, cols]))


def redundancy(weights: np.ndarray) -> float:
    normed = weights / np.maximum(np.linalg.norm(weights, axis=1, keepdims=True), 1e-12)
    similarity = normed @ normed.T
    return float(similarity[np.triu_indices(len(weights), 1)].max())


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    CACHE.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SAMPLING_SEED)

    with h5py.File(INPUT, "r") as handle:
        obs = handle["obs"]
        subtype = obs_column(obs, "subtype_new")
        patient = np.asarray([collapse_patient(x) for x in obs_column(obs, "cohort")])
        treatment = obs_column(obs, "treatment")
        cell_ids = decode(obs["_index"][:])
        genes = decode(handle["var/_index"][:])
        malignant = subtype == "cancer cells"
        n_cells = int(malignant.sum())
        n_genes = len(genes)
        x_group = handle["X"]
        indptr = x_group["indptr"][:]

        sums = np.zeros(n_genes, dtype=np.float64)
        sumsq = np.zeros(n_genes, dtype=np.float64)
        detected = np.zeros(n_genes, dtype=np.int64)
        for start in range(0, len(subtype), CHUNK_ROWS):
            end = min(start + CHUNK_ROWS, len(subtype))
            keep = malignant[start:end]
            if not keep.any():
                continue
            chunk = read_csr_chunk(x_group, indptr, start, end, n_genes)[keep]
            sums += np.asarray(chunk.sum(axis=0)).ravel()
            sumsq += np.asarray(chunk.power(2).sum(axis=0)).ravel()
            detected += np.asarray((chunk != 0).sum(axis=0)).ravel()

        mean = sums / n_cells
        variance = np.maximum(sumsq / n_cells - mean**2, 0)
        detection = detected / n_cells
        technical_exclusion = np.asarray([
            gene.upper().startswith("MT-")
            or gene.upper().startswith("RPL")
            or gene.upper().startswith("RPS")
            for gene in genes
        ])
        eligible = (
            (detection >= MIN_DETECTION)
            & (detection <= MAX_DETECTION)
            & (variance > 0)
            & ~technical_exclusion
        )
        hvg_score = standardized_variance(mean, variance, eligible)
        selected_idx = np.argsort(hvg_score)[-N_HVG:]
        selected_idx.sort()

        groups: dict[tuple[str, str], list[int]] = defaultdict(list)
        for index in np.flatnonzero(malignant):
            groups[(patient[index], treatment[index])].append(int(index))
        sampled = []
        for key in sorted(groups):
            choices = np.asarray(groups[key])
            take = min(MAX_CELLS_PER_PATIENT_TIME, len(choices))
            sampled.extend(rng.choice(choices, size=take, replace=False).tolist())
        sampled = np.asarray(sorted(sampled), dtype=np.int64)
        sampled_set = set(sampled.tolist())

        pieces = []
        sampled_ids = []
        for start in range(0, len(subtype), CHUNK_ROWS):
            end = min(start + CHUNK_ROWS, len(subtype))
            local = np.asarray([index - start for index in range(start, end) if index in sampled_set], dtype=int)
            if not len(local):
                continue
            chunk = read_csr_chunk(x_group, indptr, start, end, n_genes)
            pieces.append(chunk[local][:, selected_idx])
            sampled_ids.extend(cell_ids[start:end][local].tolist())

    x = vstack(pieces, format="csr").astype(np.float32)
    gene_scale = np.sqrt(np.maximum(variance[selected_idx], 1e-8)).astype(np.float32)
    x = x.multiply(1.0 / gene_scale).tocsr()
    save_npz(CACHE, x)

    selection = {
        "response_variables_accessed": False,
        "source_matrix": "X (nonnegative normalized/log-transformed expression)",
        "malignant_annotation": "subtype_new == cancer cells",
        "malignant_cells": n_cells,
        "detection_fraction_range": [MIN_DETECTION, MAX_DETECTION],
        "excluded_prefixes": ["MT-", "RPL", "RPS"],
        "cell_cycle_policy": "retained; annotate and model explicitly after program freeze",
        "stress_policy": "retained; annotate and model explicitly after program freeze",
        "selection_method": "top standardized log-variance within 20 equal-frequency mean-expression bins",
        "eligible_genes_before_top_n": int(eligible.sum()),
        "selected_genes": int(len(selected_idx)),
        "selected_gene_names": genes[selected_idx].tolist(),
        "balanced_sampling": "up to 100 cells per patient x treatment group; no response labels",
        "sampling_seed": SAMPLING_SEED,
        "sampled_cells": int(x.shape[0]),
        "sampled_cell_ids_sha256_pending": True,
        "gene_scaling": "divide each selected gene by malignant-cell standard deviation; no centering",
    }
    (OUTDIR / "gene_selection.json").write_text(json.dumps(selection, indent=2) + "\n")
    (OUTDIR / "sampled_cell_ids.json").write_text(json.dumps(sampled_ids, indent=2) + "\n")

    metrics = []
    weights: dict[tuple[int, int], np.ndarray] = {}
    for rank in RANKS:
        for seed in SEEDS:
            model = NMF(
                n_components=rank,
                init="random",
                random_state=seed,
                solver="cd",
                beta_loss="frobenius",
                max_iter=MAX_ITER,
                tol=1e-4,
            )
            model.fit_transform(x)
            weights[(rank, seed)] = model.components_.astype(np.float32)
            metrics.append({
                "rank": rank,
                "seed": seed,
                "reconstruction_error": float(model.reconstruction_err_),
                "iterations": int(model.n_iter_),
                "max_within_model_component_cosine": redundancy(model.components_),
            })

    for rank in RANKS:
        reference = weights[(rank, SEEDS[0])]
        for item in metrics:
            if item["rank"] == rank:
                item["stability_to_seed_11"] = matched_stability(reference, weights[(rank, item["seed"])])

    summary = []
    for rank in RANKS:
        chosen = [item for item in metrics if item["rank"] == rank]
        summary.append({
            "rank": rank,
            "mean_reconstruction_error": float(np.mean([x["reconstruction_error"] for x in chosen])),
            "sd_reconstruction_error": float(np.std([x["reconstruction_error"] for x in chosen], ddof=1)),
            "mean_stability_to_seed_11": float(np.mean([x["stability_to_seed_11"] for x in chosen])),
            "min_stability_to_seed_11": float(np.min([x["stability_to_seed_11"] for x in chosen])),
            "mean_max_redundancy": float(np.mean([x["max_within_model_component_cosine"] for x in chosen])),
        })
    model_info = {
        "response_variables_accessed": False,
        "ranks": RANKS,
        "seeds": SEEDS,
        "max_iter": MAX_ITER,
        "metrics": metrics,
        "rank_summary": summary,
    }
    (OUTDIR / "rank_seed_metrics.json").write_text(json.dumps(model_info, indent=2) + "\n")
    np.savez_compressed(
        OUTDIR / "candidate_program_weights.npz",
        genes=genes[selected_idx],
        **{f"rank_{rank}_seed_{seed}": value for (rank, seed), value in weights.items()},
    )

    ranks = np.asarray(RANKS)
    errors = np.asarray([x["mean_reconstruction_error"] for x in summary])
    stability = np.asarray([x["mean_stability_to_seed_11"] for x in summary])
    redundancy_values = np.asarray([x["mean_max_redundancy"] for x in summary])
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.6))
    axes[0].plot(ranks, errors, marker="o"); axes[0].set_title("Reconstruction error")
    axes[1].plot(ranks, stability, marker="o"); axes[1].set_title("Matched seed stability")
    axes[2].plot(ranks, redundancy_values, marker="o"); axes[2].set_title("Maximum redundancy")
    for ax in axes:
        ax.set_xlabel("NMF rank"); ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGDIR / "nmf_rank_stability.png", dpi=180)
    fig.savefig(FIGDIR / "nmf_rank_stability.pdf")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
