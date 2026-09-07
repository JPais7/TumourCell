#!/usr/bin/env python3
"""Freeze the response-blind rank-9, seed-11 malignant program definition."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix


ROOT = Path(__file__).resolve().parents[1]
DISCOVERY = ROOT / "results/phase1/blind_discovery"
WEIGHTS_INPUT = DISCOVERY / "candidate_program_weights.npz"
SELECTION_INPUT = DISCOVERY / "gene_selection.json"
MATRIX_INPUT = ROOT / "data/processed/GSE246613/blind_nmf_input.npz"
H5AD_INPUT = ROOT / "data/raw/GSE246613/GSE246613_PembroRT_non_immune_cells.h5ad"
FROZEN_NPZ = DISCOVERY / "program_definitions_v1.npz"
FROZEN_JSON = DISCOVERY / "program_definitions_v1.json"
RANK = 9
SEED = 11
TOP_N = 50

LABELS = [
    "housekeeping-cytoskeleton",
    "immediate-early-stress",
    "proliferation-cell-cycle",
    "immune-like-CD45",
    "androgen-secretory",
    "apocrine-lipid-metabolic",
    "basal-mesenchymal-ECM",
    "antigen-presentation-interferon",
    "luminal-secretory",
]

GENE_SETS = {
    "cell_cycle": {"MKI67", "TOP2A", "UBE2C", "BIRC5", "TYMS", "TK1", "TUBA1B", "TUBB", "STMN1", "CENPF", "CCNB1", "CCNB2", "CDK1"},
    "stress": {"FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "DDIT3", "HSPA1A", "HSPA1B", "HSP90AA1", "IER2", "NR4A1"},
    "immune": {"PTPRC", "CD37", "CD52", "CD7", "IL7R", "RHOH", "CST7", "LAPTM5", "RAC2", "SRGN"},
    "housekeeping": {"ACTG1", "TPT1", "PPIA", "PFN1", "H3F3A", "H3F3B", "TMSB4X", "CFL1", "NME2", "SERF2"},
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode(values) -> np.ndarray:
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values])


def obs_column(obs: h5py.Group, name: str) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = decode(obj["categories"][:])
        return np.asarray([categories[code] if code >= 0 else "" for code in obj["codes"][:]])
    return decode(obj[:])


def selected_gene_scale(selected_genes: np.ndarray) -> np.ndarray:
    with h5py.File(H5AD_INPUT, "r") as handle:
        all_genes = decode(handle["var/_index"][:])
        lookup = {gene: i for i, gene in enumerate(all_genes)}
        columns = np.asarray([lookup[gene] for gene in selected_genes])
        malignant = obs_column(handle["obs"], "subtype_new") == "cancer cells"
        matrix = handle["X"]
        indptr = matrix["indptr"][:]
        sums = np.zeros(len(columns), dtype=np.float64)
        sumsq = np.zeros(len(columns), dtype=np.float64)
        n = int(malignant.sum())
        for start in range(0, len(malignant), 2048):
            end = min(start + 2048, len(malignant))
            keep = malignant[start:end]
            if not keep.any():
                continue
            left, right = int(indptr[start]), int(indptr[end])
            chunk = csr_matrix(
                (matrix["data"][left:right], matrix["indices"][left:right], indptr[start:end + 1] - left),
                shape=(end - start, len(all_genes)),
            )[keep][:, columns]
            sums += np.asarray(chunk.sum(axis=0)).ravel()
            sumsq += np.asarray(chunk.power(2).sum(axis=0)).ravel()
    variance = np.maximum(sumsq / n - (sums / n) ** 2, 1e-8)
    return np.sqrt(variance)


def main() -> None:
    candidates = np.load(WEIGHTS_INPUT)
    genes = candidates["genes"].astype(str)
    weights = candidates[f"rank_{RANK}_seed_{SEED}"].astype(np.float64)
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), 1e-12)
    selection = json.loads(SELECTION_INPUT.read_text())
    gene_scale = selected_gene_scale(genes)

    np.savez_compressed(
        FROZEN_NPZ,
        genes=genes,
        weights=weights,
        gene_scale=gene_scale,
        labels=np.asarray(LABELS),
        rank=np.asarray(RANK),
        seed=np.asarray(SEED),
    )

    programs = []
    for index, (label, row) in enumerate(zip(LABELS, weights), start=1):
        order = np.argsort(row)[::-1]
        top = genes[order[:TOP_N]].tolist()
        overlap = {name: sorted(set(top) & members) for name, members in GENE_SETS.items()}
        programs.append({
            "program_id": f"P{index}",
            "label": label,
            "top_genes": top,
            "top_gene_weights": [float(row[i]) for i in order[:TOP_N]],
            "top_50_weight_fraction": float(row[order[:TOP_N]].sum()),
            "technical_overlap_top_50": overlap,
            "interpretation_basis": "gene weights only; response labels unopened",
        })

    manifest = {
        "version": "v1",
        "frozen": True,
        "frozen_before_response_analysis": True,
        "response_variables_accessed_during_discovery": False,
        "rank": RANK,
        "seed": SEED,
        "rank_choice": "highest minimum matched-seed stability (0.913) with near-maximal mean stability (0.960); rank 8 had a non-converged seed",
        "normalization_discovery": "H5AD X; nonnegative normalized/log-transformed cell expression; selected genes divided by malignant-cell SD without centering",
        "scoring_method": "weighted mean of pseudobulk log1p(CPM) divided by the frozen malignant-cell gene SD over 2000 genes; weights sum to one per program",
        "program_thresholds": "none; continuous scores only",
        "selected_genes": int(len(genes)),
        "gene_selection": {key: value for key, value in selection.items() if key != "selected_gene_names"},
        "programs": programs,
        "candidate_weights_sha256": sha256(WEIGHTS_INPUT),
        "frozen_npz_sha256": sha256(FROZEN_NPZ),
    }
    FROZEN_JSON.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"rank": RANK, "seed": SEED, "npz_sha256": manifest["frozen_npz_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
