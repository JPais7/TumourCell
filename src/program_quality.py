#!/usr/bin/env python3
"""Response-blind quality assessment of the frozen malignant programs."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.sparse import load_npz


ROOT = Path(__file__).resolve().parents[1]
FROZEN = ROOT / "results/phase1/blind_discovery/program_definitions_v1.npz"
INPUT = ROOT / "data/processed/GSE246613/blind_nmf_input.npz"
OUTPUT = ROOT / "results/phase1/blind_discovery/program_quality_v1.json"
TOP_COHERENCE = 30


def main() -> None:
    frozen = np.load(FROZEN)
    genes = frozen["genes"].astype(str)
    weights = frozen["weights"]
    labels = frozen["labels"].astype(str)
    x = load_npz(INPUT)
    normalized = weights / np.maximum(np.linalg.norm(weights, axis=1, keepdims=True), 1e-12)
    redundancy = normalized @ normalized.T
    cell_total = np.asarray(x.sum(axis=1)).ravel()
    records = []
    for index, (label, row) in enumerate(zip(labels, weights)):
        order = np.argsort(row)[::-1]
        top = order[:TOP_COHERENCE]
        values = x[:, top].toarray()
        corr = np.corrcoef(values, rowvar=False)
        upper = corr[np.triu_indices(TOP_COHERENCE, 1)]
        program_score = np.asarray(x @ row).ravel()
        technical_corr = float(np.corrcoef(program_score, cell_total)[0, 1])
        records.append({
            "program_id": f"P{index + 1}",
            "label": label,
            "top_30_genes": genes[top].tolist(),
            "mean_top_30_gene_pearson": float(np.nanmean(upper)),
            "median_top_30_gene_pearson": float(np.nanmedian(upper)),
            "fraction_positive_top_30_correlations": float(np.mean(upper > 0)),
            "correlation_with_selected_gene_total": technical_corr,
            "maximum_cosine_to_other_program": float(np.max(np.delete(redundancy[index], index))),
            "closest_program": f"P{np.argmax(np.where(np.arange(len(weights)) == index, -np.inf, redundancy[index])) + 1}",
        })
    payload = {
        "response_variables_accessed": False,
        "coherence_definition": "pairwise Pearson correlation among the 30 highest-weight genes in the balanced blind cell sample",
        "technical_proxy": "correlation of weighted program score with total selected-gene expression per cell",
        "programs": records,
        "program_cosine_similarity": redundancy.tolist(),
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps(records, indent=2))


if __name__ == "__main__":
    main()
