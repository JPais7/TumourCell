#!/usr/bin/env python3
"""1500/2500-feature sensitivity against frozen Phase-3 programs."""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import phase3_discover_states as discovery  # noqa: E402

COHORTS = ["ARTEMIS", "WU_GSE176078", "KARAAYVAZ_GSE118389"]
OUT = ROOT / "results/phase3b/feature_sensitivity.csv"


def select_features(counts, genes, n_features):
    lib = np.asarray(counts.sum(axis=1)).ravel()
    x = counts.multiply(10_000 / np.maximum(lib, 1)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    mean = np.asarray(x.mean(axis=0)).ravel()
    var = np.maximum(np.asarray(x.power(2).mean(axis=0)).ravel() - mean**2, 0)
    dispersion = var / np.maximum(mean, 1e-8)
    eligible = discovery.feature_mask(genes) & (detected >= max(20, int(np.ceil(.01 * x.shape[0]))))
    idx = np.flatnonzero(eligible)
    idx = idx[np.argsort(-dispersion[idx], kind="stable")[:n_features]]
    return x[:, idx].toarray().astype(np.float32), genes[idx]


def main():
    rows = []
    for name in COHORTS:
        print(name, flush=True)
        if name == "KARAAYVAZ_GSE118389": genes, counts, _, _ = discovery.karaayvaz_cohort()
        else: genes, counts, _, _ = discovery.h5ad_cohort(name)
        frozen = np.load(ROOT / f"results/phase3/cohort_programs/{name}_programs.npz")
        meta = json.loads((ROOT / f"results/phase3/cohort_programs/{name}_programs.json").read_text())
        rank = int(meta["selected_rank"])
        for nf in (1500, 2500):
            x, genes_new = select_features(counts, genes, nf)
            weights, _, _ = discovery.fit_rank(x, rank, 20260906 + nf)
            common = sorted(set(frozen["genes"]) & set(genes_new))
            fi = {g:i for i,g in enumerate(frozen["genes"])}; ni = {g:i for i,g in enumerate(genes_new)}
            corr = discovery.corr_rows(frozen["weights"][:, [fi[g] for g in common]],
                                       weights[:, [ni[g] for g in common]])
            r, q = linear_sum_assignment(-corr)
            recovered = np.zeros(len(frozen["weights"])); recovered[r] = corr[r, q]
            for i, value in enumerate(recovered):
                rows.append({"cohort": name, "candidate": f"{name}_C{i+1:02d}",
                             "n_features": nf, "common_genes": len(common),
                             "spearman_recovery": float(value), "passes_0_60": bool(value >= .60)})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


if __name__ == "__main__": main()
