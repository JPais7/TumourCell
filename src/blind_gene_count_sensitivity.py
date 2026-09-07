#!/usr/bin/env python3
"""Response-blind NMF sensitivity to selecting 1500 versus 2500 genes."""

from __future__ import annotations

import json
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import vstack
from sklearn.decomposition import NMF

from blind_program_discovery import (
    CHUNK_ROWS,
    INPUT,
    MAX_DETECTION,
    MEAN_BINS,
    MIN_DETECTION,
    decode,
    obs_column,
    read_csr_chunk,
    standardized_variance,
)


ROOT = Path(__file__).resolve().parents[1]
SAMPLED_IDS = ROOT / "results/phase1/blind_discovery/sampled_cell_ids.json"
OUTPUT_NPZ = ROOT / "results/phase1/blind_discovery/gene_count_sensitivity_weights.npz"
OUTPUT_JSON = ROOT / "results/phase1/blind_discovery/gene_count_sensitivity_metrics.json"
GENE_COUNTS = [1500, 2500]
RANK = 9
SEED = 11


def main() -> None:
    wanted_ids = set(json.loads(SAMPLED_IDS.read_text()))
    with h5py.File(INPUT, "r") as handle:
        obs = handle["obs"]
        malignant = obs_column(obs, "subtype_new") == "cancer cells"
        cell_ids = decode(obs["_index"][:])
        genes = decode(handle["var/_index"][:])
        matrix = handle["X"]
        indptr = matrix["indptr"][:]
        n_cells, n_genes = int(malignant.sum()), len(genes)
        sums = np.zeros(n_genes); sumsq = np.zeros(n_genes); detected = np.zeros(n_genes, dtype=np.int64)
        for start in range(0, len(malignant), CHUNK_ROWS):
            end = min(start + CHUNK_ROWS, len(malignant))
            keep = malignant[start:end]
            if not keep.any():
                continue
            chunk = read_csr_chunk(matrix, indptr, start, end, n_genes)[keep]
            sums += np.asarray(chunk.sum(axis=0)).ravel()
            sumsq += np.asarray(chunk.power(2).sum(axis=0)).ravel()
            detected += np.asarray((chunk != 0).sum(axis=0)).ravel()
        mean = sums / n_cells
        variance = np.maximum(sumsq / n_cells - mean**2, 0)
        detection = detected / n_cells
        excluded = np.asarray([g.upper().startswith(("MT-", "RPL", "RPS")) for g in genes])
        eligible = (detection >= MIN_DETECTION) & (detection <= MAX_DETECTION) & (variance > 0) & ~excluded
        hvg_score = standardized_variance(mean, variance, eligible)

        selections = {}
        for count in GENE_COUNTS:
            idx = np.argsort(hvg_score)[-count:]; idx.sort(); selections[count] = idx
        union = np.unique(np.concatenate(list(selections.values())))
        pieces = []
        for start in range(0, len(malignant), CHUNK_ROWS):
            end = min(start + CHUNK_ROWS, len(malignant))
            local = np.asarray([i for i, cid in enumerate(cell_ids[start:end]) if cid in wanted_ids], dtype=int)
            if len(local):
                pieces.append(read_csr_chunk(matrix, indptr, start, end, n_genes)[local][:, union])
    x_union = vstack(pieces, format="csr").astype(np.float32)
    union_lookup = {column: i for i, column in enumerate(union)}
    arrays = {}
    metrics = []
    for count, selected in selections.items():
        local = np.asarray([union_lookup[i] for i in selected])
        scale = np.sqrt(np.maximum(variance[selected], 1e-8)).astype(np.float32)
        x = x_union[:, local].multiply(1.0 / scale).tocsr()
        model = NMF(n_components=RANK, init="random", random_state=SEED, solver="cd", max_iter=500, tol=1e-4)
        model.fit_transform(x)
        arrays[f"genes_{count}"] = genes[selected]
        arrays[f"gene_scale_{count}"] = scale
        arrays[f"weights_{count}"] = model.components_.astype(np.float32)
        metrics.append({"selected_genes": count, "rank": RANK, "seed": SEED, "reconstruction_error": float(model.reconstruction_err_), "iterations": int(model.n_iter_)})
    np.savez_compressed(OUTPUT_NPZ, **arrays)
    OUTPUT_JSON.write_text(json.dumps({"response_variables_accessed": False, "models": metrics}, indent=2) + "\n")
    print(json.dumps(metrics, indent=2))


if __name__ == "__main__":
    main()
