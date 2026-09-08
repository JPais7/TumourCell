#!/usr/bin/env python3
"""Project frozen Phase-2/3 programs into HTAN TNBC malignant pseudobulks."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix

from audit_htan_htapp_mbc import read_column


def aggregate_counts(h5, selected, groups, n_groups, n_genes):
    matrix = h5["raw/X"]
    indptr = matrix["indptr"][:]
    data, indices = matrix["data"], matrix["indices"]
    counts = np.zeros((n_groups, n_genes), dtype=np.float64)
    for start in range(0, len(selected), 2048):
        end = min(start + 2048, len(selected))
        chosen = selected[start:end]
        if not chosen.any():
            continue
        left, right = int(indptr[start]), int(indptr[end])
        chunk = csr_matrix(
            (data[left:right], indices[left:right], indptr[start:end + 1] - left),
            shape=(end - start, n_genes),
        )
        local = groups[start:end]
        for index in np.unique(local[chosen]):
            counts[index] += np.asarray(chunk[local == index].sum(axis=0)).ravel()
    return counts


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--p8", type=Path, required=True)
    parser.add_argument("--csv-out", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()

    atlas = json.loads(args.atlas.read_text())
    p8 = json.loads(args.p8.read_text())
    with h5py.File(args.h5ad, "r") as h5:
        obs = h5["obs"]
        donor = read_column(obs, "donor_id")
        sample = read_column(obs, "sampleid")
        receptor = read_column(obs, "receptors_biopsy")
        compartment = read_column(obs, "compartments")
        cell_type = read_column(obs, "cell_type")
        site = read_column(obs, "site_biopsy")
        selected = ((receptor == "ER-/PR-/HER2-") &
                    (compartment == "Malignant") &
                    (cell_type == "malignant cell"))
        donors = np.asarray(sorted(set(donor[selected])))
        donor_index = {name: index for index, name in enumerate(donors)}
        groups = np.full(len(donor), -1, dtype=np.int16)
        for name, index in donor_index.items():
            groups[selected & (donor == name)] = index
        genes = read_column(h5["raw/var"], "feature_name")
        counts = aggregate_counts(h5, selected, groups, len(donors), len(genes))

    library = counts.sum(axis=1)
    gene_index = {gene: i for i, gene in enumerate(genes)}
    log_cpm = np.log1p(counts / library[:, None] * 1_000_000.0)
    rows = []
    state_scores = {}
    for state in atlas["states"]:
        present = [g for g in state["top_genes"] if g in gene_index]
        score = log_cpm[:, [gene_index[g] for g in present]].mean(axis=1)
        state_scores[state["state_id"]] = score

    present_mask = np.asarray([g in gene_index for g in p8["genes"]])
    source_indices = np.flatnonzero(present_mask)
    target_indices = np.asarray([gene_index[p8["genes"][i]] for i in source_indices])
    weights = np.asarray(p8["weights"])[source_indices]
    weights /= weights.sum()
    scales = np.asarray(p8["gene_scale"])[source_indices]
    p8_score = (log_cpm[:, target_indices] / scales) @ weights

    for i, name in enumerate(donors):
        mask = selected & (donor == name)
        row = {
            "donor_id": str(name),
            "sampleid": "|".join(sorted(set(map(str, sample[mask])))),
            "site_biopsy": "|".join(sorted(set(map(str, site[mask])))),
            "n_malignant_tnbc": int(mask.sum()),
            "library_size": int(library[i]),
            "P8": float(p8_score[i]),
        }
        row.update({state: float(values[i]) for state, values in state_scores.items()})
        rows.append(row)

    args.csv_out.parent.mkdir(parents=True, exist_ok=True)
    with args.csv_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)

    rng = np.random.default_rng(20260908)
    summary = {}
    metrics = {**state_scores, "P8": p8_score}
    for name, values in metrics.items():
        boot = np.asarray([rng.choice(values, len(values), replace=True).mean() for _ in range(10000)])
        summary[name] = {
            "mean": float(values.mean()), "sd_between_donors": float(values.std(ddof=1)),
            "range": [float(values.min()), float(values.max())],
            "bootstrap_patient_mean_ci95": [float(x) for x in np.quantile(boot, [.025, .975])],
        }
    payload = {
        "analysis": "frozen program projection into strict TNBC malignant donor pseudobulks",
        "no_refitting": True, "thresholds_optimized": False,
        "eligibility": "ER-/PR-/HER2-, Malignant, malignant cell",
        "n_donors": len(donors), "n_malignant_cells": int(selected.sum()),
        "state_score": "unweighted mean log1p(CPM) over each frozen 100-gene state set",
        "p8_score": p8["score"], "p8_present_genes": int(present_mask.sum()),
        "summary": summary,
        "interpretation_limit": "descriptive recurrence/domain-shift projection; no clinical endpoint in this object",
    }
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
