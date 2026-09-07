#!/usr/bin/env python3
"""Match outcome-blind cohort programs and freeze the construction atlas."""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.sparse import csr_matrix
from scipy.stats import rankdata

import phase3_discover_states as discovery


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/phase3"
COHORTS = ["GSE246613", "ARTEMIS", "WU_GSE176078", "KARAAYVAZ_GSE118389"]
N_NULL = 10_000
SEED = 20260906


def load_programs():
    result = {}
    x = np.load(ROOT / "results/phase1/blind_discovery/program_definitions_v1.npz")
    result["GSE246613"] = {"genes": x["genes"], "weights": x["weights"],
                            "ids": [f"GSE246613_P{i+1:02d}" for i in range(len(x["weights"]))],
                            "stability": 1.0, "cells": 6328, "platform": "10x"}
    for name in COHORTS[1:]:
        x = np.load(OUT / f"cohort_programs/{name}_programs.npz")
        meta = json.loads((OUT / f"cohort_programs/{name}_programs.json").read_text())
        selected = next(v for v in meta["screening"] if v["rank"] == meta["selected_rank"])
        result[name] = {"genes": x["genes"], "weights": x["weights"],
                        "ids": [f"{name}_C{i+1:02d}" for i in range(len(x["weights"]))],
                        "stability": selected["stability_median"], "cells": meta["metadata"]["cells_sampled"],
                        "platform": "Fluidigm/full-length" if "KARAAYVAZ" in name else "10x"}
    return result


def expression_means(programs):
    means = {}
    sparse = np.load(ROOT / "data/processed/GSE246613/blind_nmf_input.npz")
    x = csr_matrix((sparse["data"], sparse["indices"], sparse["indptr"]), shape=tuple(sparse["shape"]))
    means["GSE246613"] = dict(zip(programs["GSE246613"]["genes"], np.asarray(x.mean(axis=0)).ravel()))
    for name in COHORTS[1:]:
        if name == "KARAAYVAZ_GSE118389": genes, counts, _, _ = discovery.karaayvaz_cohort()
        else: genes, counts, _, _ = discovery.h5ad_cohort(name)
        x, chosen, _ = discovery.normalize_and_features(counts, genes)
        assert np.array_equal(chosen, programs[name]["genes"])
        means[name] = dict(zip(chosen, x.mean(axis=0)))
    return means


def spearman_matrix(a, b):
    ar = np.apply_along_axis(rankdata, 1, a); br = np.apply_along_axis(rankdata, 1, b)
    ar -= ar.mean(1, keepdims=True); br -= br.mean(1, keepdims=True)
    ar /= np.maximum(np.linalg.norm(ar, axis=1, keepdims=True), 1e-12)
    br /= np.maximum(np.linalg.norm(br, axis=1, keepdims=True), 1e-12)
    return ar @ br.T


def pair_metrics(a, b):
    common = sorted(set(a["genes"]) & set(b["genes"]))
    ia = {g:i for i,g in enumerate(a["genes"])}; ib = {g:i for i,g in enumerate(b["genes"])}
    aa = a["weights"][:, [ia[g] for g in common]]; bb = b["weights"][:, [ib[g] for g in common]]
    corr = spearman_matrix(aa, bb)
    jac = np.zeros_like(corr)
    overlap = np.zeros_like(corr, dtype=int)
    topa = [set(a["genes"][np.argsort(-w)[:50]]) for w in a["weights"]]
    topb = [set(b["genes"][np.argsort(-w)[:50]]) for w in b["weights"]]
    for i, j in itertools.product(range(len(topa)), range(len(topb))):
        overlap[i,j] = len(topa[i] & topb[j]); jac[i,j] = overlap[i,j] / len(topa[i] | topb[j])
    return common, aa, bb, corr, jac, overlap, topa, topb


def null_threshold(a, b, common, aa, bb, topa, topb, mean_a, mean_b):
    rng = np.random.default_rng(SEED + len(common) + len(a["weights"]) * 100 + len(b["weights"]))
    propensity = rankdata([np.log1p(mean_a.get(g, 0)) + np.log1p(mean_b.get(g, 0)) for g in common]) / len(common)
    bins = np.minimum((propensity * 10).astype(int), 9)
    ar = np.apply_along_axis(rankdata, 1, aa); br = np.apply_along_axis(rankdata, 1, bb)
    ar -= ar.mean(1, keepdims=True); br -= br.mean(1, keepdims=True)
    ar /= np.maximum(np.linalg.norm(ar, axis=1, keepdims=True), 1e-12)
    br /= np.maximum(np.linalg.norm(br, axis=1, keepdims=True), 1e-12)
    corr_null, jac_null = np.empty(N_NULL), np.empty(N_NULL)
    common_arr = np.asarray(common)
    for k in range(N_NULL):
        i = rng.integers(len(ar)); j = rng.integers(len(br)); perm = np.arange(len(common))
        for q in range(10):
            idx = np.flatnonzero(bins == q); perm[idx] = rng.permutation(idx)
        v = br[j, perm]; v = (v-v.mean()) / max(np.linalg.norm(v-v.mean()), 1e-12)
        corr_null[k] = ar[i] @ v
        # Apply the same expression-stratified gene-label permutation to the B top set.
        mapping = {common_arr[t]: common_arr[perm[t]] for t in range(len(common))}
        pb = {mapping.get(g, g) for g in topb[j]}
        jac_null[k] = len(topa[i] & pb) / len(topa[i] | pb)
    return float(np.quantile(corr_null, .95)), float(np.quantile(jac_null, .95))


class UnionFind:
    def __init__(self, nodes): self.parent={x:x for x in nodes}; self.cohorts={x:{x.split("_C")[0].split("_P")[0]} for x in nodes}
    def find(self,x):
        while self.parent[x]!=x: self.parent[x]=self.parent[self.parent[x]];x=self.parent[x]
        return x
    def union_safe(self,a,b,cohort_of):
        ra,rb=self.find(a),self.find(b)
        if ra==rb:return True
        ca={cohort_of[x] for x in self.parent if self.find(x)==ra};cb={cohort_of[x] for x in self.parent if self.find(x)==rb}
        if ca&cb:return False
        self.parent[rb]=ra;return True


def main():
    programs = load_programs(); means = expression_means(programs)
    node_cohort = {pid:c for c in COHORTS for pid in programs[c]["ids"]}
    pair_results, eligible_edges = [], []
    for ca, cb in itertools.combinations(COHORTS, 2):
        a,b=programs[ca],programs[cb]
        common,aa,bb,corr,jac,ov,topa,topb=pair_metrics(a,b)
        c95,j95=null_threshold(a,b,common,aa,bb,topa,topb,means[ca],means[cb])
        valid=(corr>=.30)&(jac>=.15)&(corr>c95)&(jac>j95)
        score=corr+jac
        r,q=linear_sum_assignment(-score)
        for i,j in zip(r,q):
            record={"cohort_a":ca,"candidate_a":a["ids"][i],"cohort_b":cb,"candidate_b":b["ids"][j],
                    "common_genes":len(common),"spearman":float(corr[i,j]),"jaccard_top50":float(jac[i,j]),
                    "top50_overlap":int(ov[i,j]),"null_spearman_p95":c95,"null_jaccard_p95":j95,"passes":bool(valid[i,j])}
            pair_results.append(record)
            if valid[i,j]: eligible_edges.append((float(score[i,j]),a["ids"][i],b["ids"][j],record))
    nodes=[p for c in COHORTS for p in programs[c]["ids"]]; uf=UnionFind(nodes); kept=[]
    for score,a,b,r in sorted(eligible_edges,reverse=True):
        if uf.union_safe(a,b,node_cohort): kept.append(r)
    groups={}
    for n in nodes: groups.setdefault(uf.find(n),[]).append(n)
    recurrent=[v for v in groups.values() if len(v)>=2]
    recurrent.sort(key=lambda x:(-len(x),sorted(x)))
    states=[]
    for si,members in enumerate(recurrent,1):
        sid=f"State_{si:02d}"; union=sorted(set().union(*(set(programs[node_cohort[m]]["genes"]) for m in members)))
        scores=np.zeros(len(union)); ui={g:i for i,g in enumerate(union)}
        stabs=[]
        for m in members:
            c=node_cohort[m]; k=programs[c]["ids"].index(m); w=programs[c]["weights"][k]
            ranks=rankdata(w)/len(w)
            for g,v in zip(programs[c]["genes"],ranks):scores[ui[g]]+=v
            stabs.append(programs[c]["stability"])
        scores/=len(members); order=np.argsort(-scores); top=[union[i] for i in order[:100]]
        cohorts=sorted({node_cohort[m] for m in members}); platforms={programs[c]["platform"] for c in cohorts}
        robust=len(cohorts)>=3 and len(platforms)>=2 and np.median(stabs)>=.70
        states.append({"state_id":sid,"members":members,"cohorts":cohorts,"n_cohorts":len(cohorts),
                       "top_genes":top,"median_stability":float(np.median(stabs)),"robust_construction":bool(robust)})
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/"construction_matching.json").write_text(json.dumps({"n_null":N_NULL,"seed":SEED,"pair_assignments":pair_results,"kept_edges":kept,"states":states},indent=2)+"\n")
    (OUT/"state_atlas_v1.json").write_text(json.dumps({
        "version": "1.0-construction-frozen", "outcome_blind": True,
        "matching_thresholds": {"spearman": 0.30, "jaccard_top50": 0.15,
                                "empirical_null_percentile": 95, "null_permutations": N_NULL},
        "states": states,
    }, indent=2) + "\n")
    with (ROOT/"results/state_reproducibility_matrix.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["state_id",*COHORTS,"reserved_PAL_GSE161529"])
        for s in states:w.writerow([s["state_id"],*["PRESENTE" if c in s["cohorts"] else "INDETERMINADO" for c in COHORTS],"RESERVADO"])
    with (ROOT/"results/state_metadata.csv").open("w",newline="") as f:
        w=csv.writer(f);w.writerow(["state_id","n_construction_cohorts","construction_cohorts","members","top_20_genes","median_stability","robust_construction","validation_status"])
        for s in states:w.writerow([s["state_id"],s["n_cohorts"],";".join(s["cohorts"]),";".join(s["members"]),";".join(s["top_genes"][:20]),f'{s["median_stability"]:.6f}',s["robust_construction"],"RESERVADO"])
    fig,ax=plt.subplots(figsize=(11,6));positions={}; colors=plt.cm.Set2(np.linspace(0,1,len(COHORTS)))
    for ci,c in enumerate(COHORTS):
        ids=programs[c]["ids"]
        for ni,n in enumerate(ids): positions[n]=(ci,ni-(len(ids)-1)/2);ax.scatter(ci,positions[n][1],s=80,color=colors[ci],zorder=3);ax.text(ci+.04,positions[n][1],n.split(c)[-1].strip('_'),fontsize=7,va='center')
    for e in kept:
        a,b=e["candidate_a"],e["candidate_b"];x1,y1=positions[a];x2,y2=positions[b]
        ax.plot([x1,x2],[y1,y2],color='0.35',alpha=.7,lw=1+4*e["jaccard_top50"],zorder=1)
    ax.set_xticks(range(len(COHORTS)),COHORTS,rotation=10);ax.set_ylabel("Candidato intra-coorte (ordem arbitrária)");ax.set_title("Rede congelada de correspondências entre programas malignos");ax.spines[['top','right']].set_visible(False);fig.tight_layout()
    (ROOT/"figures").mkdir(exist_ok=True);fig.savefig(ROOT/"figures/state_network.pdf");fig.savefig(ROOT/"figures/state_network.png",dpi=180);plt.close(fig)
    print(json.dumps({"states":len(states),"robust":sum(s['robust_construction'] for s in states),"kept_edges":len(kept)},indent=2))


if __name__=="__main__":main()
