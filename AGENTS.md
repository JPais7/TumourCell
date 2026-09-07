# Project instructions

## Scientific priority

Optimize for falsifiability, patient-level independence, external validation, interpretability, and biological importance. Do not treat the digital twin as the default deliverable.

## Evidence language

Label claims as fact, observation, inference, hypothesis, model assumption, or model prediction. Never infer causality from differential expression, clustering, trajectory inference, velocity, ligand–receptor scores, attention, or feature importance.

## Validation

- Never use random cell-level splitting as primary validation.
- Freeze program definitions, thresholds, covariates, and endpoints before external testing.
- Keep discovery, validation, and falsification cohorts separate.
- Cells are repeated measurements within biological samples, not independent patients.

## Data and software

- Preserve raw inputs as immutable and checksum all downloads.
- Maintain patient/sample/time/treatment manifests.
- Record software versions, seeds, configuration, exclusions, and provenance.
- Do not download a new dataset until its scientific role is recorded in `dataset_utility_matrix.csv`.

## Stop conditions

Apply `kill_criteria.md`. Negative results must be retained and reported. Do not advance to virtual perturbation unless Level 2 external prediction succeeds.
