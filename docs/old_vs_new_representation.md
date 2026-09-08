# Old versus new representation

The historical BioKey scorer and the new versioned encoder bridge are run independently over the same patient/timepoint pseudobulks using the frozen Phase 2/3 gene definitions. The audit compares scores, included patients and downstream statistical interpretation. The machine-readable comparison is `results/phase4/old_vs_new_representation.json`.

The old pipeline remains present and authoritative for the historical results until equivalence passes. No historical output is deleted or overwritten by the comparison.

The BioKey comparison passed: all 19 patients and 37 patient/timepoint rows were identical. Maximum absolute differences were below `4.5e-15` for every feature, including P8. Patient inclusion, effects and statistical conclusions are therefore unchanged.
