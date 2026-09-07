#!/usr/bin/env python3
"""Validate the frozen construction states in the reserved Pal TNBC cohort."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

import phase3_match_states as matching
from phase3_discover_pal_holdout import load_pal
import phase3_discover_states as discovery


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase3"
PAL = "PAL_GSE161529"


def main():
    frozen = json.loads((OUT / "state_atlas_v1.json").read_text())
    programs = matching.load_programs()
    px = np.load(OUT / "cohort_programs/PAL_GSE161529_programs.npz")
    pmeta = json.loads((OUT / "cohort_programs/PAL_GSE161529_programs.json").read_text())
    selected = next(x for x in pmeta["screening"] if x["rank"] == pmeta["selected_rank"])
    programs[PAL] = {"genes": px["genes"], "weights": px["weights"],
                     "ids": [f"{PAL}_C{i+1:02d}" for i in range(len(px["weights"]))],
                     "stability": selected["stability_median"], "cells": pmeta["metadata"]["cells_sampled"], "platform": "10x"}
    means = matching.expression_means(programs)
    genes, counts, _, _ = load_pal()
    x, chosen, _ = discovery.normalize_and_features(counts, genes)
    assert np.array_equal(chosen, programs[PAL]["genes"])
    means[PAL] = dict(zip(chosen, x.mean(axis=0)))

    assignments = []
    for cohort in matching.COHORTS:
        a,b=programs[cohort],programs[PAL]
        common,aa,bb,corr,jac,ov,topa,topb=matching.pair_metrics(a,b)
        c95,j95=matching.null_threshold(a,b,common,aa,bb,topa,topb,means[cohort],means[PAL])
        valid=(corr>=.30)&(jac>=.15)&(corr>c95)&(jac>j95)
        r,q=linear_sum_assignment(-(corr+jac))
        for i,j in zip(r,q):
            assignments.append({"construction_cohort":cohort,"construction_candidate":a["ids"][i],
                                "pal_candidate":b["ids"][j],"spearman":float(corr[i,j]),
                                "jaccard_top50":float(jac[i,j]),"top50_overlap":int(ov[i,j]),
                                "null_spearman_p95":c95,"null_jaccard_p95":j95,"passes":bool(valid[i,j])})
    state_results=[]
    for state in frozen["states"]:
        hits={}
        for item in assignments:
            if item["passes"] and item["construction_candidate"] in state["members"]:
                hits.setdefault(item["pal_candidate"],[]).append(item)
        best=max(hits.values(),key=len) if hits else []
        present=len(best)>=2
        state_results.append({"state_id":state["state_id"],"status":"PRESENTE" if present else "INDETERMINADO",
                              "pal_candidate":best[0]["pal_candidate"] if best else None,
                              "independent_member_matches":len(best),"matches":best,
                              "reason":"matched at least two frozen construction members" if present else
                                       "fewer than two member matches; absence not called without bootstrap power"})
    result={"construction_atlas_sha256":matching.discovery.sha256(OUT/"state_atlas_v1.json"),
            "holdout_input_sha256":pmeta["metadata"]["input_sha256"],"pal_selected_rank":pmeta["selected_rank"],
            "pal_stability_median":selected["stability_median"],"assignments":assignments,"states":state_results}
    (OUT/"pal_holdout_validation.json").write_text(json.dumps(result,indent=2)+"\n")

    cohorts=[*matching.COHORTS,PAL]
    by_state={x["state_id"]:x for x in state_results}
    with (ROOT/"results/state_reproducibility_matrix.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["state_id",*cohorts])
        for s in frozen["states"]:
            w.writerow([s["state_id"],*["PRESENTE" if c in s["cohorts"] else "INDETERMINADO" for c in matching.COHORTS],by_state[s["state_id"]]["status"]])
    with (ROOT/"results/state_metadata.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["state_id","n_construction_cohorts","construction_cohorts","members","top_20_genes","median_stability","robust_construction","pal_validation","atlas_confidence"])
        for s in frozen["states"]:
            val=by_state[s["state_id"]]["status"]
            confidence="ALTA" if s["robust_construction"] and val=="PRESENTE" else ("MODERADA" if val=="PRESENTE" else "BAIXA")
            w.writerow([s["state_id"],s["n_cohorts"],";".join(s["cohorts"]),";".join(s["members"]),";".join(s["top_genes"][:20]),f'{s["median_stability"]:.6f}',s["robust_construction"],val,confidence])
    print(json.dumps({x["state_id"]:x["status"] for x in state_results},indent=2))


if __name__=="__main__":main()
