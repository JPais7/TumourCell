# Phase 4.9.8 — Minimum dataset for clone–state longitudinal dynamics

## Search strategy

The existing breast-cancer catalogue was re-screened across GEO, SRA/ENA and EGA, prioritising paired RNA/genotype, repeated patients, treatment phases, public barcodes and reproducible code. Same patient, sample, tumour, cluster or clone was never treated as same cell or lineage.

## Candidate landscape and observability

SRP114962 is the closest candidate (Tier B): repeated treatment-associated samples with scRNA and scDNA, but no reproducible RNA-barcode↔DNA-clone table. NEOLETEXE has longitudinal treatment design but controlled metadata and no clone mapping. GSE176078 and GSE246613 provide state/RNA context without clone-resolved time. GSE205472 is a negative control without defensible malignant identity.

## Identifiability

Selection and plasticity are **NO** at the present evidence level. Population state frequencies can be described, but neither `P(state | clone, patient, timepoint)` nor a within-clone longitudinal comparison is identifiable. No direct lineage transition or causal claim is made. The highest-value missing observation is a public, cell-level RNA-barcode to validated clone/genotype mapping in repeated samples.

## Gate and safeguards

The proposed `PHASE_5A_POPULATION_DYNAMICS_ELIGIBLE` gate remains closed: E4+ malignant identity, representation transport, external prediction, treatment-conditioned transition, persistence-baseline improvement, independent validation and uncertainty modelling are not all present. Simulation and clinical recommendations remain disabled; no outcome or treatment labels define identity or features.

## Negative controls

GSE205472 (no malignant genomic evidence), GSE176078 (CNV/labels without cell-genome linkage), and GSE246613 (cross-sectional state context) are attractive but mechanistically insufficient.
