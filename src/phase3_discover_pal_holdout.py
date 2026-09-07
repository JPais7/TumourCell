#!/usr/bin/env python3
"""Discover candidates in the reserved Pal TNBC holdout without changing the atlas."""

from pathlib import Path

import h5py
import numpy as np

import phase3_discover_states as base


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE161529/breast_epithelial_integrated_holdout.h5ad"
PAL_TNBC_IDS = {"0106", "0114", "0126", "0135", "0131", "0177", "0554", "4031"}


def load_pal():
    with h5py.File(INPUT, "r") as f:
        obs = f["obs"]
        patient = base.obs_col(obs, "donor_id")
        malignant = np.asarray([
            x.startswith("pal_Patient ") and x.rsplit(" ", 1)[-1] in PAL_TNBC_IDS for x in patient
        ])
        raw_var = f["raw/var"]
        genes = base.decode(raw_var["gene_symbols"][:]) if "gene_symbols" in raw_var else base.obs_col(raw_var, "feature_name")
        candidate_rows = np.flatnonzero(malignant)
        candidate = base.read_sparse_rows(f["raw/X"], candidate_rows, len(genes))
        local_qc = ((np.asarray(candidate.sum(axis=1)).ravel() >= 500) &
                    (np.asarray((candidate > 0).sum(axis=1)).ravel() >= 300))
        rng = np.random.default_rng(base.MASTER_SEED)
        local_keep = []
        candidate_patient = patient[candidate_rows]
        for name in sorted(set(candidate_patient[local_qc])):
            idx = np.flatnonzero(local_qc & (candidate_patient == name))
            if len(idx) > base.MAX_CELLS_PER_PATIENT:
                idx = np.sort(rng.choice(idx, base.MAX_CELLS_PER_PATIENT, replace=False))
            local_keep.extend(idx.tolist())
        local_keep = np.asarray(sorted(local_keep), dtype=int)
        rows = candidate_rows[local_keep]
        matrix = candidate[local_keep]
    return genes, matrix, patient[rows], {
        "input_sha256": base.sha256(INPUT), "cells_total": int(len(patient)),
        "malignant_before_qc": int(malignant.sum()), "malignant_after_qc": int(local_qc.sum()),
        "cells_sampled": int(len(rows)), "patients_sampled": int(len(set(patient[rows]))),
        "cell_counts_by_patient": {str(p): int((patient[rows] == p).sum()) for p in sorted(set(patient[rows]))},
        "holdout_rule": "eight Pal TNBC/BRCA1-TNBC donor IDs fixed from original GEO metadata",
    }


def main():
    genes, counts, _, meta = load_pal()
    x, selected_genes, feature_meta = base.normalize_and_features(counts, genes)
    meta.update(feature_meta)
    print("PAL_GSE161529", x.shape, meta, flush=True)
    result = base.discover("PAL_GSE161529", x, selected_genes, meta)
    print("PAL_GSE161529 selected rank", result["selected_rank"], flush=True)


if __name__ == "__main__":
    main()
