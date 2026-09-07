#!/usr/bin/env python3
"""Extract sample-level metadata from an H5AD without reading expression values."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path

import h5py


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE246613/GSE246613_PembroRT_non_immune_cells.h5ad"
OUTPUT = ROOT / "data/manifests/GSE246613_nonimmune_clinical_manifest.csv"


def decode(values):
    return [x.decode("utf-8", errors="replace") if isinstance(x, bytes) else str(x) for x in values]


def obs_column(obs, name):
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = decode(obj["categories"][:])
        return [categories[int(code)] if int(code) >= 0 else "" for code in obj["codes"][:]]
    return decode(obj[:])


def unique_or_error(values, batch, field):
    found = sorted(set(values))
    if len(found) != 1:
        raise ValueError(f"{batch}: expected one {field}, found {found}")
    return found[0]


def main():
    with h5py.File(INPUT, "r") as handle:
        obs = handle["obs"]
        columns = {
            name: obs_column(obs, name)
            for name in ["batch", "cohort", "treatment", "response_group", "celltype_new", "subtype_new"]
        }

    by_batch = defaultdict(lambda: defaultdict(list))
    for i, batch in enumerate(columns["batch"]):
        for name, values in columns.items():
            if name != "batch":
                by_batch[batch][name].append(values[i])

    rows = []
    for batch, values in sorted(by_batch.items()):
        celltypes = Counter(values["celltype_new"])
        subtypes = Counter(values["subtype_new"])
        rows.append({
            "study_id": "GSE246613",
            "patient_or_tumour_id": unique_or_error(values["cohort"], batch, "cohort"),
            "batch_id": batch,
            "treatment": unique_or_error(values["treatment"], batch, "treatment"),
            "response_group": unique_or_error(values["response_group"], batch, "response_group"),
            "nonimmune_cells": len(values["cohort"]),
            "cancer_cells": subtypes.get("cancer cells", 0),
            "fibroblasts": celltypes.get("Fibroblasts", 0),
            "endothelial": celltypes.get("Endothelial", 0),
            "pvl": celltypes.get("PVL", 0),
            "epithelial_total": celltypes.get("Epithelial", 0),
            "source_file": INPUT.name,
            "source_sha256_gz": "bc727229b3401d5c8f7704945354cd36c81d1da103e290cf747c6c26da660f9a",
        })

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} biopsy/batch rows -> {OUTPUT}")


if __name__ == "__main__":
    main()
