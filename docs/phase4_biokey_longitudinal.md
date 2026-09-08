# Phase 4 — BioKey longitudinal validation

## Scope

The authenticated BioKey processed matrices contain 226,635 cells. The analysis uses only TNBC cells annotated `Cancer_cell`, aggregates raw counts for every patient/timepoint, and projects the already frozen State_01–State_04 and P8 definitions without refitting or threshold optimisation. Cells are never treated as independent replicates.

There are 33,872 TNBC malignant cells from 19 patients. At the primary eligibility threshold (at least 50 malignant cells in both Pre and On samples), cohort 1 contributes 11 paired patients and cohort 2 contributes 3. Sensitivity thresholds of 20 and 100 cells yield 18 and 13 paired patients respectively across both cohorts.

## Primary result

In cohort 1, none of the five frozen scores shows a patient-level paired change distinguishable from zero. Mean On-minus-Pre changes are State_01 −0.006, State_02 +0.031, State_03 −0.036, State_04 +0.128 and P8 −0.082; every bootstrap 95% interval crosses zero and every two-sided Wilcoxon p-value is at least 0.24.

In cohort 2, only three patients meet the primary threshold. Point estimates are State_01 +0.049, State_02 +0.281, State_03 +0.147, State_04 −0.013 and P8 +0.041. All bootstrap intervals cross zero and the sample is too small for a stable subgroup conclusion.

T-cell expansion comparisons are exploratory. `expansion` is an immune pharmacodynamic label, not pathological response or survival. No expansion/non-expansion contrast provides robust evidence after respecting patient-level replication.

## Interpretation

BioKey supplies genuine longitudinal human single-cell evidence, but it does not validate a consistent anti-PD1-induced shift in P8 or any frozen malignant state. This is useful negative/falsification evidence: it argues against treating the current programs as universal treatment-transition coordinates. It also does not rescue the digital-twin gate, because BioKey lacks a clinical outcome suitable for patient-level prospective prediction and cohort 2 has only three primary-eligible pairs.

Consequently, the project remains at an observational molecular-state-card stage. Virtual treatment recommendations remain disabled until a frozen predictor achieves externally validated patient-level discrimination with its uncertainty interval above the prespecified null threshold.

## Reproducibility

- Analysis: `src/analyze_biokey_longitudinal.py`
- Patient/timepoint scores: `results/phase4/biokey_patient_timepoint_scores.csv`
- Statistical results: `results/phase4/biokey_longitudinal_validation.json`
- Access audit: `results/phase4/biokey_access_audit.json`
- Download checksums: `data/manifests/biokey_downloads.sha256`
