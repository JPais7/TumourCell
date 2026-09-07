# Validation strategy

1. Freeze a versioned analysis specification before opening outcome labels in held-out data.
2. Preserve study-specific normalization; transfer gene-set scores rather than integrated coordinates.
3. Aggregate inference to patient/sample or fit patient-aware hierarchical models.
4. Report effect sizes and uncertainty, not cell-level p-values.
5. Require direction agreement and a minimum prespecified patient-level effect in at least one independent treatment cohort.
6. Test platform, site, response definition, and cell-depth sensitivity.
7. Treat missing genes and inaccessible raw layers as explicit transportability failures.

Primary validation is study-level external validation. Leave-one-patient-out analysis inside discovery is only a robustness check.
