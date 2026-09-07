# Architecture options

## Option 1 — Program-based discrete transition accounting (selected)

Represent each sample by patient-level distributions of interpretable malignant program usage. Compare observed timepoints within patient and, where credible, within CNA clone. Advantages: auditable, compatible with sparse data, easy to transport. Limitation: does not reconstruct individual cell paths.

## Option 2 — Markov/state-space model

Potentially useful with several repeated timepoints and adequate patients per regimen. Current clinical data are too sparse for unconstrained transition matrices; parameters would be prior-dominated. Reject for v0.1.

## Option 3 — RNA velocity/CellRank

Can provide exploratory local directionality if raw spliced/unspliced counts exist and kinetics pass diagnostics. Treatment, patient, and tumor heterogeneity violate simple assumptions. Keep as non-primary sensitivity analysis.

## Option 4 — Neural ODE or hybrid mechanistic/ML twin

Too flexible and non-identifiable for current observations. Reject until repeated patient measurements and external prediction succeed.

## v0.1 data flow

Raw immutable inputs → metadata audit → study-specific QC → malignant/CNA evidence → locked functional programs → patient/sample summaries → discrete change and nondetection model → external validation → falsification/sensitivity report.
