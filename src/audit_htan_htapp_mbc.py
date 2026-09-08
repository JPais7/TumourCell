#!/usr/bin/env python3
"""Audit HTAN/HTAPP MBC metadata and frozen-program gene coverage.

This script deliberately does not score, fit, cluster, or alter any frozen state.
It reads the CELLxGENE H5AD directly with h5py so the audit has minimal
dependencies and writes only cohort/sample summaries.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import h5py
import numpy as np


def decode(value):
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def read_column(obs, name):
    obj = obs[name]
    if isinstance(obj, h5py.Group) and obj.attrs.get("encoding-type") == "categorical":
        categories = np.asarray([decode(x) for x in obj["categories"][()]])
        codes = obj["codes"][()]
        out = np.full(codes.shape, "<NA>", dtype=object)
        valid = codes >= 0
        out[valid] = categories[codes[valid]]
        return out
    values = obj[()]
    if values.dtype.kind in {"O", "S", "U"}:
        return np.asarray([decode(x) for x in values], dtype=object)
    return values


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", type=Path, required=True)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--p8", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    parser.add_argument("--manifest-out", type=Path, required=True)
    args = parser.parse_args()

    atlas = json.loads(args.atlas.read_text())
    p8 = json.loads(args.p8.read_text())

    with h5py.File(args.h5ad, "r") as h5:
        obs = h5["obs"]
        columns = {
            name: read_column(obs, name)
            for name in (
                "donor_id", "sampleid", "receptors_biopsy", "site_biopsy",
                "compartments", "cell_type", "author_cell_type", "assay",
                "suspension_type", "tissue", "Phase",
            )
        }
        genes = set(map(str, read_column(h5["var"], "feature_name")))
        n_cells, n_genes = map(int, h5["X"].attrs["shape"])

    tnbc = columns["receptors_biopsy"] == "ER-/PR-/HER2-"
    malignant = (columns["compartments"] == "Malignant") & (
        columns["cell_type"] == "malignant cell"
    )

    sample_rows = []
    for sample in sorted(set(columns["sampleid"])):
        mask = columns["sampleid"] == sample
        values = {}
        for field in (
            "donor_id", "receptors_biopsy", "site_biopsy", "assay",
            "suspension_type", "tissue", "Phase",
        ):
            unique = sorted(set(map(str, columns[field][mask])))
            values[field] = "|".join(unique)
        sample_rows.append({
            "sampleid": sample,
            **values,
            "n_cells": int(mask.sum()),
            "n_malignant": int((mask & malignant).sum()),
            "tnbc_eligible": bool((mask & tnbc).any()),
        })

    state_coverage = []
    for state in atlas["states"]:
        state_genes = state["top_genes"]
        present = [gene for gene in state_genes if gene in genes]
        state_coverage.append({
            "state_id": state["state_id"],
            "n_frozen_top_genes": len(state_genes),
            "n_present": len(present),
            "fraction_present": len(present) / len(state_genes),
            "missing": [gene for gene in state_genes if gene not in genes],
        })

    p8_genes = p8["genes"]
    p8_present = [gene for gene in p8_genes if gene in genes]
    result = {
        "audit_version": "1.0",
        "input": str(args.h5ad),
        "input_bytes": args.h5ad.stat().st_size,
        "input_sha256": sha256(args.h5ad),
        "shape": {"cells": n_cells, "genes": n_genes},
        "cohort": {
            "donors": len(set(columns["donor_id"])),
            "samples": len(set(columns["sampleid"])),
            "malignant_cells": int(malignant.sum()),
            "receptor_counts_cells": dict(Counter(map(str, columns["receptors_biopsy"]))),
            "site_counts_cells": dict(Counter(map(str, columns["site_biopsy"]))),
        },
        "tnbc": {
            "definition": "receptors_biopsy == ER-/PR-/HER2-",
            "cells": int(tnbc.sum()),
            "malignant_cells": int((tnbc & malignant).sum()),
            "donors": len(set(columns["donor_id"][tnbc])),
            "samples": len(set(columns["sampleid"][tnbc])),
            "donor_malignant_counts": dict(sorted(Counter(
                map(str, columns["donor_id"][tnbc & malignant])
            ).items())),
            "site_counts_malignant": dict(Counter(
                map(str, columns["site_biopsy"][tnbc & malignant])
            )),
        },
        "frozen_state_top_gene_coverage": state_coverage,
        "p8_coverage": {
            "n_frozen_genes": len(p8_genes),
            "n_present": len(p8_present),
            "fraction_present": len(p8_present) / len(p8_genes),
            "missing": [gene for gene in p8_genes if gene not in genes],
        },
        "no_model_fitting_or_scoring_performed": True,
    }

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(result, indent=2) + "\n")
    with args.manifest_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sample_rows[0]))
        writer.writeheader()
        writer.writerows(sample_rows)


if __name__ == "__main__":
    main()
