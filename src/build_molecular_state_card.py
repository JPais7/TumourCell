#!/usr/bin/env python3
"""Build a non-causal molecular state card from frozen-signature scores."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata


SIGNATURES = ["State_01", "State_02", "State_03", "State_04", "P8"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scores", type=Path, required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--dataset", choices=["GSE260693", "NeoTRIP"], required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    with args.scores.open() as handle:
        rows = list(csv.DictReader(handle))
    matches = [row for row in rows if row["sample_id"] == args.sample_id]
    if len(matches) != 1:
        raise SystemExit(f"Expected one sample named {args.sample_id!r}; found {len(matches)}")
    target = matches[0]
    reference = [row for row in rows if row["timepoint"] == target["timepoint"]]
    percentiles = {}
    for signature in SIGNATURES:
        values = np.asarray([float(row[signature]) for row in reference])
        target_index = next(i for i, row in enumerate(reference) if row["sample_id"] == args.sample_id)
        percentiles[signature] = float((rankdata(values, method="average")[target_index] - .5) / len(values) * 100)

    if args.dataset == "GSE260693":
        missing = ["cnv", "mutation", "clone", "spatial", "ctdna", "imaging"]
        treatment = target.get("nac_simple") or None
        domain = "bulk RNA-seq; scores mix malignant and microenvironmental expression"
    else:
        missing = ["clinical", "treatment", "outcome", "cnv", "mutation", "clone", "spatial", "ctdna", "imaging"]
        treatment = None
        domain = "bulk RNA-seq; clinical arm unavailable and scores mix malignant and microenvironmental expression"
    card = {
        "patient_id": target["patient_id"], "sample_id": target["sample_id"],
        "timepoint": target["timepoint"], "assay": "bulk_rna",
        "treatment_context": treatment,
        "program_scores": {signature: float(target[signature]) for signature in SIGNATURES},
        "reference_percentiles": percentiles,
        "uncertainty": {
            "measurement": f"single tissue sample; percentile referenced to {len(reference)} {args.dataset} samples at the same timepoint",
            "domain_shift": domain,
            "biological": "selection, plasticity and sampling-region effects cannot be separated",
        },
        "missing_modalities": missing,
        "allowed_uses": ["observational_state_description", "cohort_level_research"],
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(card, indent=2) + "\n")


if __name__ == "__main__":
    main()
