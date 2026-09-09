# Phase 4 implementation audit

The Phase 4 engineering pass adds a first probabilistic research model without activating a Digital Twin. The model is `p(z_next | z_current, treatment, delta_t)` with Gaussian residual covariance, explicit treatment and time interval, patient-level duplicate rejection, model hash and leave-one-patient-out validation metrics. `predict_mean`, `predict_distribution`, `sample_next` and `log_likelihood` are separate APIs; sampling is research-only and is not a treatment recommendation.

Uncertainty is represented as measurement uncertainty, sampling uncertainty, coverage penalty, extrapolation flag, missing required modalities and an overall quality state. The aggregated `uncertainty_score` is a heuristic quality score, not a calibrated confidence interval or posterior uncertainty. Missing required modality is fail-closed and cannot receive `PASS`.

The spatial layer now supports validated coordinates, k-nearest-neighbour and radius graphs, isolated-cell handling, local state summaries, neighbourhood composition and a descriptive spatial-confounding result. It does not infer causal niches or implement a spatial simulator. The observation model is an interface only.

The mechanism layer distinguishes observed state shift, composition shift, selection, plasticity, environmental shift and sampling/technical confounding. Clone evidence alone is insufficient: assignment quality, temporal comparability, abundance and sampling/QC are required. Otherwise the result is `NOT_IDENTIFIABLE` or `INDETERMINATE`.

Validation: `pytest -q` reports 30 passed and 0 failed. Historical Phase 1–3 definitions and outputs were not modified. The BioKey old/new equivalence remains passed at numerical precision. G1–G10 remain conservative; G2, G4, G5, G6, G8 and G10 are not sufficient to unlock simulation. `simulation_allowed()` returns `False` by design.
