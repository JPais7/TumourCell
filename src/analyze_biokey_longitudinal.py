#!/usr/bin/env python3
"""Patient-level longitudinal projection of frozen malignant programs in BioKey."""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd
import rdata
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu, wilcoxon


def ci_mean(values, rng, n=10000):
    values = np.asarray(values, float)
    boot = np.asarray([rng.choice(values, len(values), replace=True).mean() for _ in range(n)])
    return [float(x) for x in np.quantile(boot, [.025, .975])]


def read_sparse_rds(path):
    obj = rdata.read_rds(path)
    dim = tuple(int(x) for x in obj.Dim)
    matrix = csc_matrix((obj.x, obj.i, obj.p), shape=dim)
    return matrix, np.asarray(obj.Dimnames[0], str), np.asarray(obj.Dimnames[1], str)


def process_cohort(cohort, root, wanted_genes):
    metadata = pd.read_csv(root / f"BIOKEY_metaData_cohort{cohort}_web.csv")
    matrix, genes, cells = read_sparse_rds(root / f"counts_cells_cohort{cohort}.rds")
    if matrix.shape != (len(genes), len(cells)):
        raise ValueError(f"Cohort {cohort}: inconsistent matrix dimensions")
    meta = metadata.set_index("Cell", drop=False)
    if not np.array_equal(cells, metadata["Cell"].astype(str).to_numpy()):
        missing = set(cells) - set(meta.index)
        if missing:
            raise ValueError(f"Cohort {cohort}: {len(missing)} matrix cells absent from metadata")
        metadata = meta.loc[cells].reset_index(drop=True)

    selected = (metadata.BC_type.eq("TNBC") & metadata.cellType.eq("Cancer_cell")).to_numpy()
    gene_index = {g: i for i, g in enumerate(genes)}
    present = [g for g in wanted_genes if g in gene_index]
    small = matrix[[gene_index[g] for g in present], :]
    rows, counts = [], []
    for (patient, timepoint), frame in metadata[selected].groupby(["patient_id", "timepoint"]):
        idx = frame.index.to_numpy()
        rows.append({
            "cohort": int(cohort), "patient_id": patient, "timepoint": timepoint,
            "expansion": frame.expansion.dropna().astype(str).iloc[0] if frame.expansion.notna().any() else None,
            "n_malignant_cells": int(len(idx)), "library_size": int(matrix[:, idx].sum()),
        })
        counts.append(np.asarray(small[:, idx].sum(axis=1)).ravel())
    del matrix, small, genes
    gc.collect()
    return pd.DataFrame(rows), np.asarray(counts), present


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data", type=Path, default=Path("data/raw/BIOKEY"))
    p.add_argument("--atlas", type=Path, default=Path("results/phase3/state_atlas_v1.json"))
    p.add_argument("--p8", type=Path, default=Path("results/phase2/p8_definition_frozen.json"))
    p.add_argument("--out", type=Path, default=Path("results/phase4"))
    args = p.parse_args()
    atlas, p8 = json.loads(args.atlas.read_text()), json.loads(args.p8.read_text())
    state_genes = {s["state_id"]: s["top_genes"] for s in atlas["states"]}
    wanted = list(dict.fromkeys(g for gs in state_genes.values() for g in gs))
    wanted += [g for g in p8["genes"] if g not in wanted]

    frames = []
    for cohort in (1, 2):
        info, counts, present = process_cohort(cohort, args.data, wanted)
        program_lib = counts.sum(axis=1)
        lib = info["library_size"].to_numpy(float)
        log_cpm = np.log1p(counts / np.maximum(lib[:, None], 1) * 1e6)
        gi = {g: i for i, g in enumerate(present)}
        for state, genes in state_genes.items():
            idx = [gi[g] for g in genes if g in gi]
            info[state] = log_cpm[:, idx].mean(axis=1)
        source = [i for i, g in enumerate(p8["genes"]) if g in gi]
        weights = np.asarray(p8["weights"], float)[source]
        weights /= weights.sum()
        scales = np.asarray(p8["gene_scale"], float)[source]
        info["P8"] = (log_cpm[:, [gi[p8["genes"][i]] for i in source]] / scales) @ weights
        info["program_library_size"] = program_lib.astype(int)
        info["p8_present_genes"] = len(source)
        frames.append(info)
        del counts, log_cpm
        gc.collect()
    scores = pd.concat(frames, ignore_index=True)
    args.out.mkdir(parents=True, exist_ok=True)
    scores.to_csv(args.out / "biokey_patient_timepoint_scores.csv", index=False)

    rng = np.random.default_rng(20260908)
    metrics = list(state_genes) + ["P8"]
    analyses = {}
    for threshold in (20, 50, 100):
        threshold_results = {}
        for cohort in (1, 2):
            wide_n = scores[scores.cohort.eq(cohort)].pivot(index="patient_id", columns="timepoint", values="n_malignant_cells")
            eligible = wide_n.index[(wide_n.get("Pre", 0) >= threshold) & (wide_n.get("On", 0) >= threshold)]
            cres = {"n_paired_patients": int(len(eligible)), "metrics": {}}
            for metric in metrics:
                wide = scores[scores.cohort.eq(cohort) & scores.patient_id.isin(eligible)].pivot(index="patient_id", columns="timepoint", values=metric)
                delta = (wide.On - wide.Pre).dropna()
                stat = wilcoxon(delta).pvalue if len(delta) >= 3 and np.any(delta != 0) else None
                item = {"mean_delta_on_minus_pre": float(delta.mean()) if len(delta) else None,
                        "bootstrap_patient_mean_delta_ci95": ci_mean(delta, rng) if len(delta) else None,
                        "wilcoxon_p_two_sided": float(stat) if stat is not None else None}
                exp = scores[scores.cohort.eq(cohort)].drop_duplicates("patient_id").set_index("patient_id").expansion
                de, dn = delta[exp.reindex(delta.index).eq("E")], delta[exp.reindex(delta.index).eq("NE")]
                item["expansion_subgroup"] = {"n_E": len(de), "n_NE": len(dn), "mean_delta_E": float(de.mean()) if len(de) else None,
                    "mean_delta_NE": float(dn.mean()) if len(dn) else None,
                    "mannwhitney_p_two_sided": float(mannwhitneyu(de, dn).pvalue) if len(de) >= 2 and len(dn) >= 2 else None}
                cres["metrics"][metric] = item
            threshold_results[f"cohort_{cohort}"] = cres
        analyses[str(threshold)] = threshold_results
    payload = {
        "analysis": "BioKey longitudinal TNBC malignant pseudobulk projection",
        "unit_of_inference": "patient", "no_refitting": True, "thresholds_optimized": False,
        "primary_minimum_cells_per_timepoint": 50, "sensitivity_thresholds": [20, 100],
        "cohort_interpretation": {"1": "treatment-naive anti-PD1 context", "2": "prior chemotherapy/neoadjuvant chemotherapy plus anti-PD1 context"},
        "n_tnbc_malignant_cells": int(scores.n_malignant_cells.sum()),
        "n_tnbc_patients": int(scores.patient_id.nunique()),
        "score_definition": "patient-timepoint summed malignant counts; log1p(CPM); frozen state means and frozen scaled P8 weighted score",
        "analyses": analyses,
        "interpretation_limit": "Expansion is an immune pharmacodynamic subgroup, not pathological response or survival; results do not authorize treatment simulation.",
    }
    (args.out / "biokey_longitudinal_validation.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
