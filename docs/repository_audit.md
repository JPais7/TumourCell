# Repository audit

Audit date: 2026-09-08. Historical Phase 1–3/3c outputs remain unchanged. `tests/` and `experiments/` existed but were empty; tests are now executable and experiment manifests have a defined writer. `notebooks/` remains empty because scripts are the source of truth. Raw/processed data remain ignored. The principal debt is migration of each historical scorer to the new encoder after locked old-vs-new comparisons; the old scripts remain available until equivalence is documented.

No available dataset currently identifies clone-resolved selection versus plasticity or validates a patient-specific counterfactual. These remain explicit stop conditions.
