#!/usr/bin/env python3
"""Project the frozen Phase-2/3 signatures into GSE260693 and NeoTRIP."""

from __future__ import annotations

import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import chi2, mannwhitneyu, norm, ttest_1samp, wilcoxon


ROOT = Path(__file__).resolve().parents[1]
ATLAS = json.loads((ROOT / "results/phase3/state_atlas_v1.json").read_text())
P8 = json.loads((ROOT / "results/phase2/p8_definition_frozen.json").read_text())
OUT = ROOT / "results/phase3d"
OUT.mkdir(parents=True, exist_ok=True)


def geo_metadata(path: Path) -> list[dict[str, str]]:
    fields: dict[str, list[str]] = {}
    repeated = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                fields["sample_title"] = [x.strip('"') for x in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_geo_accession"):
                fields["geo_accession"] = [x.strip('"') for x in next(csv.reader([line], delimiter="\t"))[1:]]
            elif line.startswith("!Sample_characteristics_ch1"):
                values = [x.strip('"') for x in next(csv.reader([line], delimiter="\t"))[1:]]
                key = values[0].split(":", 1)[0]
                repeated[key] = [x.split(":", 1)[1].strip() if ":" in x else x.strip() for x in values]
    fields.update(repeated)
    return [{key: values[i] for key, values in fields.items()} for i in range(len(fields["sample_title"]))]


def read_expression(path: Path, delimiter: str, leading_columns: int, log1p: bool) -> tuple[list[str], dict[str, np.ndarray]]:
    expression: dict[str, list[np.ndarray]] = defaultdict(list)
    wanted = set(P8["genes"])
    for state in ATLAS["states"]:
        wanted.update(state["top_genes"])
    with gzip.open(path, "rt") as handle:
        reader = csv.reader(handle, delimiter=delimiter)
        header = next(reader)
        samples = [x.strip('"') for x in header[leading_columns:]]
        for row in reader:
            if not row:
                continue
            gene = row[0].strip('"')
            if gene not in wanted:
                continue
            values = np.asarray([float(x) if x not in {"", "NA", "NaN"} else np.nan for x in row[leading_columns:]])
            if log1p:
                values = np.log1p(np.maximum(values, 0))
            expression[gene].append(values)
    return samples, {gene: np.nanmean(values, axis=0) for gene, values in expression.items()}


def score(expression: dict[str, np.ndarray]) -> tuple[dict[str, np.ndarray], dict[str, dict[str, float]]]:
    scores, coverage = {}, {}
    for state in ATLAS["states"]:
        present = [gene for gene in state["top_genes"] if gene in expression]
        scores[state["state_id"]] = np.nanmean([expression[gene] for gene in present], axis=0)
        coverage[state["state_id"]] = {"present": len(present), "total": len(state["top_genes"])}
    indices = [i for i, gene in enumerate(P8["genes"]) if gene in expression]
    weights_raw = np.asarray(P8["weights"], dtype=float)[indices]
    weights = weights_raw / weights_raw.sum()
    scales = np.asarray(P8["gene_scale"], dtype=float)[indices]
    matrix = np.asarray([expression[P8["genes"][i]] for i in indices]).T / scales
    available = np.isfinite(matrix)
    scores["P8"] = np.nansum(matrix * weights, axis=1) / (available @ weights)
    coverage["P8"] = {
        "present": len(indices), "total": len(P8["genes"]),
        "weight_present": float(weights_raw.sum()),
    }
    return scores, coverage


def bootstrap_difference(values: np.ndarray, groups: np.ndarray, seed: int) -> dict:
    a, b = values[groups], values[~groups]
    rng = np.random.default_rng(seed)
    draws = np.asarray([
        rng.choice(a, len(a), replace=True).mean() - rng.choice(b, len(b), replace=True).mean()
        for _ in range(10000)
    ])
    pooled = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    return {
        "n_group_1": len(a), "n_group_0": len(b),
        "mean_difference": float(a.mean() - b.mean()),
        "bootstrap_ci95": [float(x) for x in np.quantile(draws, [.025, .975])],
        "hedges_g": float((a.mean() - b.mean()) / pooled * (1 - 3 / (4 * (len(a) + len(b)) - 9))),
        "mann_whitney_p_two_sided": float(mannwhitneyu(a, b, alternative="two-sided").pvalue),
    }


def paired_test(delta: np.ndarray, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    draws = np.asarray([rng.choice(delta, len(delta), replace=True).mean() for _ in range(10000)])
    try:
        p_wilcoxon = float(wilcoxon(delta, alternative="two-sided").pvalue)
    except ValueError:
        p_wilcoxon = 1.0
    return {
        "n_pairs": len(delta), "mean_delta_later_minus_baseline": float(delta.mean()),
        "median_delta": float(np.median(delta)),
        "bootstrap_mean_ci95": [float(x) for x in np.quantile(draws, [.025, .975])],
        "wilcoxon_p_two_sided": p_wilcoxon,
        "paired_t_p_two_sided": float(ttest_1samp(delta, 0).pvalue),
    }


def bh(tests: dict, p_key: str) -> None:
    keys = list(tests)
    p = np.asarray([tests[key][p_key] for key in keys])
    order = np.argsort(p)
    adjusted = np.empty(len(p))
    running = 1.0
    for rank in range(len(p) - 1, -1, -1):
        idx = order[rank]
        running = min(running, p[idx] * len(p) / (rank + 1))
        adjusted[idx] = running
    for key, value in zip(keys, adjusted):
        tests[key]["bh_fdr"] = float(value)


def write_scores(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze_gse260693() -> None:
    metadata = geo_metadata(ROOT / "tmp/GSE260693/series_matrix.txt.gz")
    by_title = {row["sample_title"]: row for row in metadata}
    samples, expression = read_expression(
        ROOT / "tmp/GSE260693/GSE260693_tpm_normalized_counts.csv.gz", ",", 1, True
    )
    scores, coverage = score(expression)
    rows = []
    for i, sample in enumerate(samples):
        meta = by_title[sample]
        timepoint = sample.split("_", 1)[0]
        row = {
            "sample_id": sample, "patient_id": meta["patient id"], "timepoint": timepoint,
            "residual_tumor": meta["residual tumor"], "nac_simple": meta["nac simple"],
        }
        row.update({name: float(values[i]) for name, values in scores.items()})
        rows.append(row)
    write_scores(OUT / "gse260693_frozen_scores.csv", rows)

    baseline_idx = np.asarray([row["timepoint"] == "pre" for row in rows])
    residual = np.asarray([row["residual_tumor"].lower() == "yes" for row in rows])[baseline_idx]
    outcome_tests = {
        name: bootstrap_difference(values[baseline_idx], residual, 20260908 + i)
        for i, (name, values) in enumerate(scores.items())
    }
    bh(outcome_tests, "mann_whitney_p_two_sided")

    positions = {(row["patient_id"], row["timepoint"]): i for i, row in enumerate(rows)}
    paired_patients = sorted({p for p, t in positions if t == "pre" and (p, "post") in positions})
    paired_tests = {}
    for i, (name, values) in enumerate(scores.items()):
        delta = np.asarray([values[positions[(p, "post")]] - values[positions[(p, "pre")]] for p in paired_patients])
        paired_tests[name] = paired_test(delta, 20261008 + i)
    bh(paired_tests, "wilcoxon_p_two_sided")
    payload = {
        "analysis": "GSE260693 frozen signatures; log1p(TPM); no refitting",
        "samples": len(samples), "unique_patients": len({row["patient_id"] for row in rows}),
        "baseline": int(baseline_idx.sum()), "post_treatment": int((~baseline_idx).sum()),
        "complete_pairs": len(paired_patients), "baseline_no_residual": int((~residual).sum()),
        "baseline_residual": int(residual.sum()), "coverage": coverage,
        "outcome_contrast": "residual tumour yes minus no at baseline",
        "baseline_outcome_tests": outcome_tests, "paired_change_tests": paired_tests,
        "limitations": ["bulk tissue composition", "heterogeneous NAC regimens", "19 complete pairs", "outcome is residual tumour rather than centrally adjudicated pCR"],
    }
    (OUT / "gse260693_validation.json").write_text(json.dumps(payload, indent=2) + "\n")


def analyze_neotrip() -> None:
    samples, expression = read_expression(ROOT / "tmp/GSE319641/GSE319641_TPM.txt.gz", "\t", 3, False)
    scores, coverage = score(expression)
    rows = []
    for i, sample in enumerate(samples):
        patient, timepoint = sample.rsplit("_", 1)
        row = {"sample_id": sample, "patient_id": patient, "timepoint": timepoint}
        row.update({name: float(values[i]) for name, values in scores.items()})
        rows.append(row)
    write_scores(OUT / "neotrip_frozen_scores.csv", rows)
    positions = {(row["patient_id"], row["timepoint"]): i for i, row in enumerate(rows)}
    paired_patients = sorted({p for p, t in positions if t == "Baseline" and (p, "D1C2") in positions})
    tests = {}
    for i, (name, values) in enumerate(scores.items()):
        delta = np.asarray([values[positions[(p, "D1C2")]] - values[positions[(p, "Baseline")]] for p in paired_patients])
        tests[name] = paired_test(delta, 20261108 + i)
    bh(tests, "wilcoxon_p_two_sided")
    payload = {
        "analysis": "NeoTRIP GSE319641 frozen signature longitudinal projection; no refitting",
        "samples": len(samples), "unique_patients": len({row["patient_id"] for row in rows}),
        "baseline": sum(row["timepoint"] == "Baseline" for row in rows),
        "d1c2": sum(row["timepoint"] == "D1C2" for row in rows),
        "complete_pairs": len(paired_patients), "coverage": coverage, "paired_change_tests": tests,
        "clinical_endpoint_tested": False,
        "limitations": ["clinical pCR and treatment arm are not public", "bulk tissue composition", "ComBat-adjusted TPM supplied by authors"],
    }
    (OUT / "neotrip_longitudinal.json").write_text(json.dumps(payload, indent=2) + "\n")


def fixed_effect_meta() -> None:
    ispy = json.loads((ROOT / "results/phase3c/ispy2_tnbc_validation.json").read_text())
    gse = json.loads((OUT / "gse260693_validation.json").read_text())
    results = {}
    for signature in [state["state_id"] for state in ATLAS["states"]] + ["P8"]:
        source = ispy["unadjusted_tests"][signature]
        validation = gse["baseline_outcome_tests"][signature]
        # Orient every effect as favourable outcome (pCR/no residual) minus residual disease.
        effects = np.asarray([source["hedges_g"], -validation["hedges_g"]])
        sample_sizes = [(source["n_pcr"], source["n_residual"]),
                        (validation["n_group_0"], validation["n_group_1"])]
        variances = np.asarray([
            (n1 + n0) / (n1 * n0) + effect**2 / (2 * (n1 + n0 - 2))
            for effect, (n1, n0) in zip(effects, sample_sizes)
        ])
        weights = 1 / variances
        combined = float(np.sum(weights * effects) / weights.sum())
        se = float(np.sqrt(1 / weights.sum()))
        q = float(np.sum(weights * (effects - combined) ** 2))
        results[signature] = {
            "orientation": "pCR/no residual minus residual disease",
            "cohort_hedges_g": {"ISPY2": float(effects[0]), "GSE260693": float(effects[1])},
            "fixed_effect_hedges_g": combined, "standard_error": se,
            "ci95": [combined - 1.959964 * se, combined + 1.959964 * se],
            "p_two_sided": float(2 * norm.sf(abs(combined / se))),
            "cochran_q": q, "heterogeneity_p": float(chi2.sf(q, 1)),
            "i2_percent": float(max(0, (q - 1) / q) * 100) if q else 0.0,
        }
    bh(results, "p_two_sided")
    payload = {
        "analysis": "fixed-effect meta-analysis of two independent pretreatment TNBC cohorts",
        "cohorts": ["I-SPY2", "GSE260693"], "no_refitting": True,
        "effect_orientation": "favourable pathological outcome minus residual disease",
        "signatures": results,
        "caution": "only two cohorts; report individual effects and heterogeneity alongside pooled estimates",
    }
    (OUT / "ispy2_gse260693_meta_analysis.json").write_text(json.dumps(payload, indent=2) + "\n")


if __name__ == "__main__":
    analyze_gse260693()
    analyze_neotrip()
    fixed_effect_meta()
