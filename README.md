# Dynamic Tumour Cell

Scientific reconnaissance for a patient-informed model of breast-cancer state dynamics. The project deliberately starts with a falsifiable biological question, not a full digital twin.

## Current decision

**Proceed with a TNBC v0.2 focused on longitudinal malignant-program redistribution during pembrolizumab followed by pembrolizumab plus radiotherapy.** Do not build a patient-specific simulator or claim that selection and plasticity have been separated. Public data support treatment-associated state change, but sparse clinical sampling does not identify continuous-time transition rates.

The primary reconnaissance deliverable is [docs/reconnaissance_report.md](docs/reconnaissance_report.md). The pre-analysis contract is in [docs/phase1_analysis_specification.md](docs/phase1_analysis_specification.md), the access gate is in [docs/metadata_access_audit.md](docs/metadata_access_audit.md), dataset-level evidence is in [dataset_utility_matrix.csv](dataset_utility_matrix.csv), and the recommended portfolio is in [minimum_sufficient_portfolio.md](minimum_sufficient_portfolio.md).

## Proposed v0.1 question

Which malignant-cell programs change from baseline to pembrolizumab and from pembrolizumab to pembrolizumab plus radiotherapy, and which changes distinguish early responders, late responders and nonresponders?

## Scope

- Unit of inference: patient/sample, never cells treated as independent replicates.
- Discovery: program-level malignant states, with stress and cell cycle modeled explicitly.
- Dynamics: discrete observed transitions only (pre, on-treatment, residual), not continuous trajectories unless raw spliced/unspliced data and sampling support them.
- Validation: locked program definitions and thresholds in independent studies.
- Output: an externally tested biological prediction and a calibrated statement of what the data cannot identify.

## Repository layout

`data/`, `src/`, `notebooks/`, `experiments/`, `results/`, `figures/`, `docs/`, `tests/`, and `configs/` are reserved for Phase 1. Raw biological files are excluded from version control; manifests and provenance remain tracked.

## Status

Reconnaissance completed 2026-09-06. The discovery gate is approved for malignant cells in GSE246613. The primary depth threshold is frozen at 50 malignant cells at each of three timepoints, with mandatory sensitivity analyses. Kim/SRP114962 remains a future clone-aware extension. Machine-generated manifests are under `data/manifests/`, and the first count-level pseudobulk is under `data/processed/GSE246613/`.
