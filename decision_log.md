# Decision log

## 2026-09-06 — Do not build the digital twin

Reason: public data support discrete state comparisons, not patient-specific continuous dynamics.

## 2026-09-06 — Select TNBC conditionally

Reason: strongest combination of single-cell, treatment, response, spatial, and independent cohorts. Reconsider if raw/metadata access fails.

## 2026-09-06 — Reject broad pretreatment response classifier

Reason: Yan et al. 2026 already provides a larger, integrated single-cell/spatial ecotype analysis and response panel.

## 2026-09-06 — Select program usage over hard clusters

Reason: mixed functional programs are more interpretable and transportable; clusters are sensitive to batch and resolution.

## 2026-09-06 — Select Portfolio B

Reason: adds distinct information at each layer without the redundancy and harmonization burden of the high-confidence portfolio.

## 2026-09-06 — Gate 1 partially passes after metadata audit

Kim raw data are public but operationally unsuitable for full first-pass reprocessing (~1,386 GiB across 6,633 runs) without a processed matrix and cell dictionary. GSE246613 has 34 molecularly profiled patients and strong serial completeness, but response labels are absent from GEO SOFT. GSE205472 has only five complete pairs and remains an exploratory control.

## 2026-09-06 — GSE246613 clinical metadata recovered

The human non-immune H5AD contains response group, treatment, patient/tumour identifier and malignant annotations for 171,746 cells. After collapsing two tumours from Patient03, the cohort has 34 patients: 9 R1, 14 R2 and 11 NR. Malignant-cell depth across all three times is adequate for only 12 patients at a 50-cell threshold, so threshold sensitivity is mandatory.

## 2026-09-06 — TISCH2 is evidence, not yet a usable source

TISCH2 lists a processed BRCA_SRP114962 object with 2,472 cells from eight patients, but its public download link returned 404. The project will not claim the matrix is available until the bytes and metadata dictionary are retrieved and checksummed.

## 2026-09-06 — GSE169246 does not replace Kim for malignant dynamics

The public RNA barcode file contains 489,490 CD45+ cells across 78 samples from 22 patients. Twelve patients have paired pre/post tumour samples, evenly split between paclitaxel and paclitaxel plus atezolizumab, but the CD45 enrichment removes the malignant compartment required by the primary question. The cohort is retained only for immune-mechanism validation.

## 2026-09-06 — Pivot primary discovery to GSE246613

Kim remains inaccessible as an auditable processed object, and GSE169246 contains only CD45+ cells. The project therefore pivots from clone-aware chemotherapy dynamics to longitudinal malignant-program redistribution during pembrolizumab followed by pembrolizumab plus radiotherapy. The primary 50-cell-per-timepoint threshold is frozen before expression inspection, all depth sensitivities remain mandatory, and no claim separating selection from plasticity is allowed.

## 2026-09-06 — Freeze nine programs before response analysis

Response-blind NMF rank 9, seed 11 was selected from ranks 4–10 and five fixed seeds using reconstruction, matched-seed stability and redundancy. The frozen definition hash is c353b2b5edb47c18e960aeb10622dc673e6f6fdf24342deebd784a3982371aec. No response variable was loaded by the discovery or quality-control scripts.

## 2026-09-06 — Stop complexity escalation after Phase 1

The primary cohort contains six NR, five R2 and one R1 patient. No R2-versus-NR program contrast has a bootstrap interval excluding zero. P8 antigen-presentation/interferon retains a positive direction across ranks, seeds, gene counts, normalizations and cell-depth thresholds, but remains imprecise and can reverse when Patient45 is removed. Selection versus plasticity is not identifiable because the H5AD contains no CNV or clonal linkage. Do not build a digital twin or response predictor; proceed only to partial external validation and a CNV-feasibility audit.
