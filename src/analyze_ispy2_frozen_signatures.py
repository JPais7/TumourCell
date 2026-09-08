#!/usr/bin/env python3
"""Project frozen states/P8 into I-SPY2 pretreatment bulk expression."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu, t as student_t


def sample_metadata(paths):
    result = {}
    for path in paths:
        fields = {}
        with gzip.open(path, "rt") as handle:
            for line in handle:
                if line.startswith("!Sample_title"):
                    row = next(csv.reader([line], delimiter="\t"))[1:]
                    fields["title"] = [x.strip('"') for x in row]
                elif line.startswith("!Sample_characteristics_ch1"):
                    row = [x.strip('"') for x in next(csv.reader([line], delimiter="\t"))[1:]]
                    key = row[0].split(":", 1)[0]
                    fields[key] = [x.split(":", 1)[1].strip() for x in row]
        for i, title in enumerate(fields["title"]):
            patient = fields["patient id"][i]
            result[patient] = {k: v[i] for k, v in fields.items() if k != "title"}
            result[patient]["sample_title"] = title
            result[patient]["platform_source"] = path.name
    return result


def bootstrap_difference(values, outcome, seed=20260908):
    a, b = values[outcome == 1], values[outcome == 0]
    rng = np.random.default_rng(seed)
    boot = np.asarray([
        rng.choice(a, len(a), replace=True).mean() - rng.choice(b, len(b), replace=True).mean()
        for _ in range(10000)
    ])
    pooled = np.sqrt(((len(a)-1)*a.var(ddof=1)+(len(b)-1)*b.var(ddof=1))/(len(a)+len(b)-2))
    return {
        "n_pcr": len(a), "n_residual": len(b),
        "mean_difference_pcr_minus_residual": float(a.mean()-b.mean()),
        "bootstrap_ci95": [float(x) for x in np.quantile(boot, [.025, .975])],
        "hedges_g": float((a.mean()-b.mean())/pooled*(1-3/(4*(len(a)+len(b))-9))),
        "mann_whitney_p_two_sided": float(mannwhitneyu(a, b, alternative="two-sided").pvalue),
    }


def adjusted_pcr_effect(values, outcome, meta):
    arms = np.asarray([m["arm"] for m in meta])
    platforms = np.asarray([m["platform_source"] for m in meta])
    columns = [np.ones(len(values)), outcome.astype(float)]
    for level in sorted(set(arms))[1:]:
        columns.append((arms == level).astype(float))
    for level in sorted(set(platforms))[1:]:
        columns.append((platforms == level).astype(float))
    design = np.column_stack(columns)
    beta, _, _, _ = np.linalg.lstsq(design, values, rcond=None)
    residual = values - design @ beta
    dof = len(values) - design.shape[1]
    covariance = np.linalg.inv(design.T @ design) * (residual @ residual / dof)
    se = float(np.sqrt(covariance[1, 1]))
    statistic = float(beta[1] / se)
    critical = float(student_t.ppf(.975, dof))
    return {"coefficient_pcr_adjusted_arm_platform": float(beta[1]),
            "standard_error": se, "ci95": [float(beta[1]-critical*se), float(beta[1]+critical*se)],
            "p_two_sided": float(2*student_t.sf(abs(statistic), dof)), "degrees_freedom": dof}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--expression", type=Path, required=True)
    parser.add_argument("--series-matrix", type=Path, nargs="+", required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--p8", type=Path, required=True)
    parser.add_argument("--scores-out", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()
    atlas = json.loads(args.atlas.read_text())
    p8 = json.loads(args.p8.read_text())
    metadata = sample_metadata(args.series_matrix)
    wanted = set(p8["genes"])
    for state in atlas["states"]:
        wanted.update(state["top_genes"])
    expression = defaultdict(list)
    with gzip.open(args.expression, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        patients = next(reader)
        for row in reader:
            if row and row[0] in wanted:
                expression[row[0]].append(np.asarray([
                    float(value) if value not in {"NA", "", "NaN"} else np.nan
                    for value in row[1:]
                ]))
    # Gene-level file can contain repeated symbols; average them without choosing
    # a response-favourable probe.
    expression = {g: np.nanmean(v, axis=0) for g, v in expression.items()}
    keep = np.asarray([p in metadata for p in patients])
    patients = np.asarray(patients)[keep]
    meta = [metadata[p] for p in patients]
    tnbc = np.asarray([m["hr"] == "0" and m["her2"] == "0" for m in meta])
    outcome = np.asarray([int(m["pcr"]) for m in meta])[tnbc]
    selected_patients = patients[tnbc]

    scores = {}
    coverage = {}
    for state in atlas["states"]:
        genes = state["top_genes"]
        present = [g for g in genes if g in expression]
        coverage[state["state_id"]] = {"present": len(present), "total": len(genes)}
        scores[state["state_id"]] = np.nanmean(
            [expression[g][keep][tnbc] for g in present], axis=0
        )
    present_indices = [i for i, g in enumerate(p8["genes"]) if g in expression]
    weights = np.asarray(p8["weights"])[present_indices]
    weights /= weights.sum()
    scales = np.asarray(p8["gene_scale"])[present_indices]
    matrix = np.asarray([expression[p8["genes"][i]][keep][tnbc] for i in present_indices]).T
    scaled = matrix / scales
    available = np.isfinite(scaled)
    scores["P8"] = np.nansum(scaled * weights, axis=1) / (available @ weights)
    coverage["P8"] = {"present": len(present_indices), "total": len(p8["genes"]),
                       "weight_present": float(np.asarray(p8["weights"])[present_indices].sum())}

    rows = []
    for i, patient in enumerate(selected_patients):
        m = metadata[patient]
        row = {"patient_id": patient, "pcr": m["pcr"], "arm": m["arm"]}
        row.update({name: float(value[i]) for name, value in scores.items()})
        rows.append(row)
    args.scores_out.parent.mkdir(parents=True, exist_ok=True)
    with args.scores_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)

    tests = {name: bootstrap_difference(value, outcome, 20260908+i)
             for i, (name, value) in enumerate(scores.items())}
    selected_meta = [m for m, use in zip(meta, tnbc) if use]
    for name, value in scores.items():
        tests[name]["adjusted"] = adjusted_pcr_effect(value, outcome, selected_meta)
    payload = {
        "analysis": "frozen signatures in I-SPY2 TNBC pretreatment bulk microarray",
        "no_refitting": True, "multiple_testing_correction": "Benjamini-Hochberg across five prespecified signatures",
        "patients_expression": int(keep.sum()), "tnbc_patients": int(tnbc.sum()),
        "tnbc_pcr": int(outcome.sum()), "tnbc_residual": int((outcome == 0).sum()),
        "coverage": coverage, "unadjusted_tests": tests,
        "limitations": ["bulk tissue composition", "microarray-to-scRNA domain shift", "arm heterogeneity"],
    }
    pvals = np.asarray([tests[k]["mann_whitney_p_two_sided"] for k in tests])
    order = np.argsort(pvals); adjusted = np.empty(len(pvals)); running = 1.0
    for rank in range(len(pvals)-1, -1, -1):
        idx = order[rank]; running = min(running, pvals[idx]*len(pvals)/(rank+1)); adjusted[idx] = running
    for key, value in zip(tests, adjusted): tests[key]["bh_fdr"] = float(value)
    args.json_out.write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    main()
