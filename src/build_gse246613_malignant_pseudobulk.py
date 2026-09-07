#!/usr/bin/env python3
"""Aggregate malignant raw counts by patient and treatment without loading all cells."""

from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "data/raw/GSE246613/GSE246613_PembroRT_non_immune_cells.h5ad"
OUTPUT = ROOT / "data/processed/GSE246613/malignant_pseudobulk_counts.npz"
QC_OUTPUT = ROOT / "results/phase1/gse246613_malignant_qc.json"
TREATMENT_ORDER = ["Base", "PD1", "RTPD1"]
THRESHOLDS = [1, 10, 25, 50, 100]


def decode(values) -> np.ndarray:
    return np.asarray([
        value.decode("utf-8", errors="replace") if isinstance(value, bytes) else str(value)
        for value in values
    ])


def obs_column(obs: h5py.Group, name: str) -> np.ndarray:
    obj = obs[name]
    if isinstance(obj, h5py.Group) and "categories" in obj and "codes" in obj:
        categories = decode(obj["categories"][:])
        codes = obj["codes"][:]
        return np.asarray([categories[code] if code >= 0 else "" for code in codes])
    return decode(obj[:])


def collapse_patient(value: str) -> str:
    return "Patient03" if value in {"Patient03T1", "Patient03T2"} else value


def main() -> None:
    with h5py.File(INPUT, "r") as handle:
        obs = handle["obs"]
        subtype = obs_column(obs, "subtype_new")
        tumour = obs_column(obs, "cohort")
        patient = np.asarray([collapse_patient(value) for value in tumour])
        treatment = obs_column(obs, "treatment")
        response = obs_column(obs, "response_group")
        batches = obs_column(obs, "batch")
        genes = decode(handle["raw/var/_index"][:])

        malignant = subtype == "cancer cells"
        group_keys = sorted(
            set(zip(patient[malignant], treatment[malignant])),
            key=lambda item: (item[0], TREATMENT_ORDER.index(item[1])),
        )
        group_index = {key: i for i, key in enumerate(group_keys)}
        cell_group = np.full(len(subtype), -1, dtype=np.int16)
        for key, index in group_index.items():
            cell_group[malignant & (patient == key[0]) & (treatment == key[1])] = index

        counts = np.zeros((len(group_keys), len(genes)), dtype=np.float64)
        matrix = handle["raw/X"]
        indptr = matrix["indptr"][:]
        data_ds, indices_ds = matrix["data"], matrix["indices"]
        chunk_rows = 2048
        for start in range(0, len(subtype), chunk_rows):
            end = min(start + chunk_rows, len(subtype))
            selected = cell_group[start:end] >= 0
            if not selected.any():
                continue
            left, right = int(indptr[start]), int(indptr[end])
            chunk = csr_matrix(
                (
                    data_ds[left:right],
                    indices_ds[left:right],
                    indptr[start : end + 1] - left,
                ),
                shape=(end - start, len(genes)),
            )
            local_groups = cell_group[start:end]
            for index in np.unique(local_groups[selected]):
                counts[index] += np.asarray(chunk[local_groups == index].sum(axis=0)).ravel()

    group_patients = np.asarray([key[0] for key in group_keys])
    group_treatments = np.asarray([key[1] for key in group_keys])
    group_cells = np.asarray([(malignant & (patient == p) & (treatment == t)).sum() for p, t in group_keys])
    response_by_patient = {}
    for name in sorted(set(patient[malignant])):
        found = sorted(set(response[malignant & (patient == name)]))
        if len(found) != 1:
            raise ValueError(f"{name}: inconsistent response groups {found}")
        response_by_patient[name] = found[0]
    group_responses = np.asarray([response_by_patient[name] for name in group_patients])

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        OUTPUT,
        counts=counts,
        genes=genes,
        patient=group_patients,
        treatment=group_treatments,
        response_group=group_responses,
        malignant_cells=group_cells,
    )

    depth = defaultdict(dict)
    for (name, time), n_cells in zip(group_keys, group_cells):
        depth[name][time] = int(n_cells)
    eligible = {
        str(threshold): sorted(
            name for name, values in depth.items()
            if all(values.get(time, 0) >= threshold for time in TREATMENT_ORDER)
        )
        for threshold in THRESHOLDS
    }
    qc = {
        "source_file": INPUT.name,
        "source_sha256_gz": "bc727229b3401d5c8f7704945354cd36c81d1da103e290cf747c6c26da660f9a",
        "expression_layer": "raw/X integer counts",
        "genes": int(len(genes)),
        "malignant_cells": int(malignant.sum()),
        "patients_with_malignant_cells": int(len(depth)),
        "pseudobulk_groups": int(len(group_keys)),
        "response_counts": dict(Counter(response_by_patient.values())),
        "eligible_all_three_times": {key: len(value) for key, value in eligible.items()},
        "eligible_patient_ids": eligible,
        "patient03_policy": "Patient03T1 and Patient03T2 collapsed to Patient03",
        "batch_count_malignant": int(len(set(batches[malignant]))),
        "counts_reconciliation": int(counts.sum()),
    }
    QC_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    QC_OUTPUT.write_text(json.dumps(qc, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(qc["eligible_all_three_times"], sort_keys=True))
    print(f"{counts.shape} pseudobulk matrix -> {OUTPUT}")


if __name__ == "__main__":
    main()
