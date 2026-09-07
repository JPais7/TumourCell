#!/usr/bin/env python3
"""Outcome-blind, study-specific malignant-state discovery for Phase 3."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import h5py
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.stats import rankdata
from sklearn.decomposition import MiniBatchNMF


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase3/cohort_programs"
SEEDS = [11, 23, 47, 71, 101]
RANKS = list(range(5, 16))
MASTER_SEED = 20260906
MAX_CELLS_PER_PATIENT = 2000
N_FEATURES = 2000

INPUTS = {
    "ARTEMIS": ROOT / "data/raw/ARTEMIS/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad",
    "WU_GSE176078": ROOT / "data/raw/GSE176078/GSE176078_cellxgene.h5ad",
    "KARAAYVAZ_GSE118389": ROOT / "data/raw/GSE118389/GSE118389_counts_rsem.txt.gz",
}


def decode(values):
    return np.asarray([x.decode() if isinstance(x, bytes) else str(x) for x in values])


def obs_col(group, name):
    x = group[name]
    if isinstance(x, h5py.Group):
        cats, codes = decode(x["categories"][:]), x["codes"][:]
        return np.asarray([cats[i] if i >= 0 else "<NA>" for i in codes])
    return decode(x[:])


def sha256(path, block=8 * 1024 * 1024):
    h = hashlib.sha256()
    with path.open("rb") as f:
        while b := f.read(block): h.update(b)
    return h.hexdigest()


def balanced_sample(mask, patient, qc):
    rng = np.random.default_rng(MASTER_SEED)
    keep = []
    for p in sorted(set(patient[mask & qc])):
        idx = np.flatnonzero(mask & qc & (patient == p))
        if len(idx) > MAX_CELLS_PER_PATIENT:
            idx = np.sort(rng.choice(idx, MAX_CELLS_PER_PATIENT, replace=False))
        keep.extend(idx.tolist())
    return np.asarray(sorted(keep), dtype=np.int64)


def read_sparse_rows(group, rows, n_genes, selected_cols=None):
    """Read chosen CSR rows without materializing the full matrix."""
    indptr = group["indptr"][:]
    data, indices = group["data"], group["indices"]
    blocks = []
    for start in range(0, len(rows), 1024):
        rr = rows[start : start + 1024]
        first, last = int(rr[0]), int(rr[-1]) + 1
        a, b = int(indptr[first]), int(indptr[last])
        block = csr_matrix((data[a:b], indices[a:b], indptr[first:last+1] - a), shape=(last-first, n_genes))
        block = block[rr - first]
        if selected_cols is not None: block = block[:, selected_cols]
        blocks.append(block)
    from scipy.sparse import vstack
    return vstack(blocks, format="csr")


def h5ad_cohort(name):
    path = INPUTS[name]
    with h5py.File(path, "r") as f:
        o = f["obs"]
        patient = obs_col(o, "donor_id")
        if name == "ARTEMIS":
            malignant = obs_col(o, "author_cell_type") == "Tumor"
        else:
            malignant = (obs_col(o, "subtype") == "TNBC") & (obs_col(o, "celltype_major") == "Cancer Epithelial")
        ncount = np.asarray(o["nCount_RNA"][:])
        nfeature = np.asarray(o["nFeature_RNA"][:])
        qc = (ncount >= 500) & (nfeature >= 300)
        rows = balanced_sample(malignant, patient, qc)
        raw_var = f["raw/var"]
        genes = decode(raw_var["gene_symbols"][:]) if "gene_symbols" in raw_var else obs_col(raw_var, "feature_name")
        matrix = read_sparse_rows(f["raw/X"], rows, len(genes))
    return genes, matrix, patient[rows], {
        "input_sha256": sha256(path), "cells_total": int(len(patient)),
        "malignant_before_qc": int(malignant.sum()), "malignant_after_qc": int((malignant & qc).sum()),
        "cells_sampled": int(len(rows)), "patients_sampled": int(len(set(patient[rows]))),
        "cell_counts_by_patient": {str(p): int((patient[rows] == p).sum()) for p in sorted(set(patient[rows]))},
    }


def karaayvaz_cohort():
    path = INPUTS["KARAAYVAZ_GSE118389"]
    with gzip.open(path, "rt") as f:
        cells = f.readline().rstrip("\n").split("\t")
        genes, rows = [], []
        for line in f:
            fields = line.rstrip("\n").split("\t")
            genes.append(fields[0]); rows.append(np.asarray(fields[1:], dtype=np.float32))
    values = np.vstack(rows).T
    patient = np.asarray([cell.split("_")[0] for cell in cells])
    library = values.sum(axis=1)
    detected = (values > 0).sum(axis=1)
    qc = (library >= 500) & (detected >= 300)
    values, patient = values[qc], patient[qc]
    return np.asarray(genes), csr_matrix(values), patient, {
        "input_sha256": sha256(path), "cells_total": len(cells), "malignant_before_qc": len(cells),
        "malignant_after_qc": int(qc.sum()), "cells_sampled": int(qc.sum()),
        "patients_sampled": int(len(set(patient))),
        "cell_counts_by_patient": {str(p): int((patient == p).sum()) for p in sorted(set(patient))},
        "measurement_warning": "RSEM expected counts are fractional and not UMI counts",
    }


def feature_mask(genes):
    bad = re.compile(r"^(MT-|RPL\d|RPS\d|IG[HKL][VDJCMAGET]|HSP[A-Z0-9]|HIST[0-9]|MALAT1$)")
    cycle = {"MKI67","PCNA","TOP2A","UBE2C","CCNB1","CCNB2","CDK1","CDC20","TYMS","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7"}
    return np.asarray([not bad.match(g) and g not in cycle for g in genes])


def normalize_and_features(counts, genes):
    lib = np.asarray(counts.sum(axis=1)).ravel()
    x = counts.multiply(10_000 / np.maximum(lib, 1)[:, None]).tocsr()
    x.data = np.log1p(x.data)
    detected = np.asarray((counts > 0).sum(axis=0)).ravel()
    mean = np.asarray(x.mean(axis=0)).ravel()
    mean2 = np.asarray(x.power(2).mean(axis=0)).ravel()
    var = np.maximum(mean2 - mean**2, 0)
    dispersion = var / np.maximum(mean, 1e-8)
    eligible = feature_mask(genes) & (detected >= max(20, int(np.ceil(.01 * x.shape[0]))))
    idx = np.flatnonzero(eligible)
    idx = idx[np.argsort(-dispersion[idx], kind="stable")[:N_FEATURES]]
    return x[:, idx].toarray().astype(np.float32), genes[idx], {
        "genes_input": int(len(genes)), "genes_eligible": int(eligible.sum()), "features_selected": int(len(idx))
    }


def corr_rows(a, b):
    ar = np.apply_along_axis(rankdata, 1, a)
    br = np.apply_along_axis(rankdata, 1, b)
    ar -= ar.mean(axis=1, keepdims=True); br -= br.mean(axis=1, keepdims=True)
    ar /= np.maximum(np.linalg.norm(ar, axis=1, keepdims=True), 1e-12)
    br /= np.maximum(np.linalg.norm(br, axis=1, keepdims=True), 1e-12)
    return ar @ br.T


def fit_rank(x, rank, seed):
    # Random initialization makes the frozen seeds a genuine stability test.
    # tol=0 prevents the one-epoch parameter-tolerance stop observed with NNDSVD;
    # max_no_improvement supplies the explicit mini-batch stopping rule.
    model = MiniBatchNMF(n_components=rank, init="random", random_state=seed, max_iter=20,
                         batch_size=2048, beta_loss="frobenius", tol=0.0,
                         max_no_improvement=None, fresh_restarts=True)
    usage = model.fit_transform(x)
    weights = model.components_
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), 1e-12)
    return weights.astype(np.float32), float(model.reconstruction_err_), int(model.n_iter_)


def stability(weights):
    ref = weights[0]
    vals = []
    for other in weights[1:]:
        c = corr_rows(ref, other)
        r, q = linear_sum_assignment(-c)
        vals.extend(c[r, q].tolist())
    return float(np.median(vals)), float(np.quantile(vals, .1))


def discover(name, x, genes, meta):
    screening = []
    fits = {}
    for rank in RANKS:
        ww, errors, iters = [], [], []
        for seed in SEEDS:
            w, e, n = fit_rank(x, rank, seed)
            ww.append(w); errors.append(e); iters.append(n)
        med, q10 = stability(ww)
        fits[rank] = ww
        screening.append({"rank": rank, "stability_median": med, "stability_q10": q10,
                          "reconstruction_error_mean": float(np.mean(errors)), "iterations": iters})
    max_stability = max(v["stability_median"] for v in screening)
    candidates = [v for v in screening if v["stability_median"] >= .99 * max_stability]
    selected = min(candidates, key=lambda v: v["rank"])["rank"]
    ref = fits[selected][0]
    aligned = [ref]
    for other in fits[selected][1:]:
        c = corr_rows(ref, other); r, q = linear_sum_assignment(-c)
        order = np.empty(selected, dtype=int); order[r] = q
        aligned.append(other[order])
    consensus = np.mean(aligned, axis=0)
    consensus /= consensus.sum(axis=1, keepdims=True)
    top = []
    for i, row in enumerate(consensus):
        order = np.argsort(-row)[:100]
        top.append({"candidate": f"{name}_C{i+1:02d}", "top_genes": genes[order].tolist(),
                    "top_weights": row[order].astype(float).tolist()})
    payload = {"cohort": name, "selected_rank": selected, "rank_rule": "smallest rank within 1% of maximal median seed stability",
               "nmf": {"implementation": "sklearn MiniBatchNMF", "init": "random", "epochs": 20,
                       "batch_size": 2048, "tol": 0.0, "max_no_improvement": None},
               "seeds": SEEDS, "ranks": RANKS, "metadata": meta, "screening": screening, "programs": top}
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / f"{name}_programs.npz", genes=genes, weights=consensus)
    (OUT / f"{name}_programs.json").write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort", choices=[*INPUTS, "ALL"], default="ALL")
    args = parser.parse_args()
    names = list(INPUTS) if args.cohort == "ALL" else [args.cohort]
    for name in names:
        if name == "KARAAYVAZ_GSE118389": genes, counts, patient, meta = karaayvaz_cohort()
        else: genes, counts, patient, meta = h5ad_cohort(name)
        x, selected_genes, feature_meta = normalize_and_features(counts, genes)
        meta.update(feature_meta)
        print(name, x.shape, meta, flush=True)
        result = discover(name, x, selected_genes, meta)
        print(name, "selected rank", result["selected_rank"], flush=True)


if __name__ == "__main__": main()
