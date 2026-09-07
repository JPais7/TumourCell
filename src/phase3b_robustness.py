#!/usr/bin/env python3
"""Diagnostic robustness audit of the frozen Phase-3 construction programs.

This script never changes the atlas. It asks whether frozen within-cohort programs
can be recovered after patient bootstrap, cell subsampling, adjacent ranks and
removal of their ten highest-weight genes.
"""

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
import phase3_match_states as matching  # noqa: E402

OUT = ROOT / "results/phase3b"
COHORTS = ["ARTEMIS", "WU_GSE176078", "KARAAYVAZ_GSE118389"]
BOOTSTRAPS = 20
CELL_REPLICATES = 5
RECOVERY_THRESHOLD = 0.60
SEED = 20260906


def align_recovery(reference, candidate):
    corr = discovery.corr_rows(reference, candidate)
    r, q = linear_sum_assignment(-corr)
    # At rank k-1 at least one frozen component can be left unmatched.  This is
    # a genuine failure of recovery, represented explicitly as zero rather than
    # NaN so that summaries and JSON remain standards-compliant.
    values = np.zeros(len(reference), dtype=float)
    values[r] = corr[r, q]
    return values


def load_cohort(name):
    if name == "KARAAYVAZ_GSE118389":
        genes, counts, patient, _ = discovery.karaayvaz_cohort()
    else:
        genes, counts, patient, _ = discovery.h5ad_cohort(name)
    x, chosen, _ = discovery.normalize_and_features(counts, genes)
    frozen = np.load(ROOT / f"results/phase3/cohort_programs/{name}_programs.npz")
    if not np.array_equal(chosen, frozen["genes"]):
        raise RuntimeError(f"Feature mismatch for {name}")
    meta = json.loads((ROOT / f"results/phase3/cohort_programs/{name}_programs.json").read_text())
    return x, patient, frozen["weights"], int(meta["selected_rank"])


def fit_and_recover(x, reference, rank, seed):
    weights, _, _ = discovery.fit_rank(x, rank, seed)
    return align_recovery(reference, weights)


def audit_cohort(name):
    x, patient, reference, rank = load_cohort(name)
    rng = np.random.default_rng(SEED + sum(map(ord, name)))
    ids = [f"{name}_C{i+1:02d}" for i in range(len(reference))]
    records = []

    unique = np.asarray(sorted(set(patient)))
    for rep in range(BOOTSTRAPS):
        sampled_patients = rng.choice(unique, len(unique), replace=True)
        rows = np.concatenate([np.flatnonzero(patient == p) for p in sampled_patients])
        rec = fit_and_recover(x[rows], reference, rank, SEED + 1000 + rep)
        for cid, value in zip(ids, rec):
            records.append({"cohort": name, "candidate": cid, "test": "patient_bootstrap",
                            "level": "100% patients with replacement", "replicate": rep + 1,
                            "spearman_recovery": float(value)})

    for fraction in (0.25, 0.50, 0.75):
        for rep in range(CELL_REPLICATES):
            rows = []
            for p in unique:
                ix = np.flatnonzero(patient == p)
                n = max(1, int(np.ceil(fraction * len(ix))))
                rows.extend(rng.choice(ix, n, replace=False).tolist())
            rec = fit_and_recover(x[np.asarray(rows)], reference, rank, SEED + int(100*fraction) + rep)
            for cid, value in zip(ids, rec):
                records.append({"cohort": name, "candidate": cid, "test": "cell_subsample",
                                "level": f"{int(100*fraction)}% per patient", "replicate": rep + 1,
                                "spearman_recovery": float(value)})

    for adjacent in sorted({max(2, rank - 1), rank + 1}):
        rec = fit_and_recover(x, reference, adjacent, SEED + adjacent)
        for cid, value in zip(ids, rec):
            records.append({"cohort": name, "candidate": cid, "test": "adjacent_rank",
                            "level": f"rank={adjacent}", "replicate": 1,
                            "spearman_recovery": float(value)})
    return records


def top10_edge_audit():
    programs = matching.load_programs()
    atlas = json.loads((ROOT / "results/phase3/state_atlas_v1.json").read_text())
    cohort_of = {pid: cohort for cohort, obj in programs.items() for pid in obj["ids"]}
    rows = []
    for state in atlas["states"]:
        members = state["members"]
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a_id, b_id = members[i], members[j]
                ca, cb = cohort_of[a_id], cohort_of[b_id]
                a, b = programs[ca], programs[cb]
                ia, ib = a["ids"].index(a_id), b["ids"].index(b_id)
                common = sorted(set(a["genes"]) & set(b["genes"]))
                amap = {g:k for k,g in enumerate(a["genes"])}
                bmap = {g:k for k,g in enumerate(b["genes"])}
                remove = set(a["genes"][np.argsort(-a["weights"][ia])[:10]]) | set(
                    b["genes"][np.argsort(-b["weights"][ib])[:10]])
                kept = [g for g in common if g not in remove]
                aw = a["weights"][[ia]][:, [amap[g] for g in kept]]
                bw = b["weights"][[ib]][:, [bmap[g] for g in kept]]
                rho = float(matching.spearman_matrix(aw, bw)[0, 0])
                rows.append({"state_id": state["state_id"], "candidate_a": a_id,
                             "candidate_b": b_id, "genes_remaining": len(kept),
                             "spearman_without_union_top10": rho,
                             "passes_rho_0_30": rho >= 0.30})
    return rows


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    records = []
    for cohort in COHORTS:
        print(f"Auditing {cohort}", flush=True)
        records.extend(audit_cohort(cohort))
    with (OUT / "program_recovery.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(records[0])); w.writeheader(); w.writerows(records)
    edge_rows = top10_edge_audit()
    with (OUT / "top10_removed_edges.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(edge_rows[0])); w.writeheader(); w.writerows(edge_rows)
    summary = []
    for cohort in COHORTS:
        candidates = sorted({r["candidate"] for r in records if r["cohort"] == cohort})
        for candidate in candidates:
            for test in ("patient_bootstrap", "cell_subsample", "adjacent_rank"):
                vals = [r["spearman_recovery"] for r in records if r["candidate"] == candidate and r["test"] == test]
                summary.append({"cohort": cohort, "candidate": candidate, "test": test,
                                "n": len(vals), "median": float(np.median(vals)),
                                "q10": float(np.quantile(vals, .1)),
                                "pass_fraction": float(np.mean(np.asarray(vals) >= RECOVERY_THRESHOLD))})
    payload = {"diagnostic_only": True, "atlas_modified": False,
               "recovery_threshold": RECOVERY_THRESHOLD, "patient_bootstraps": BOOTSTRAPS,
               "cell_replicates_per_fraction": CELL_REPLICATES, "summary": summary,
               "top10_removed_edges": edge_rows}
    (OUT / "robustness_summary.json").write_text(json.dumps(payload, indent=2) + "\n")
    print(json.dumps({"records": len(records), "edge_tests": len(edge_rows)}, indent=2))


if __name__ == "__main__":
    main()
