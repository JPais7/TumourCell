# Minimum sufficient dataset portfolio

## Portfolio A — Minimal

1. Kim 2018: longitudinal, clone-aware TNBC discovery.
2. Yan/ARTEMIS 2026: large independent pretreatment response and spatial context.
3. Wu 2021: independent cross-subtype program recurrence.

Coverage is strong for malignant programs, weak for independent post-treatment replication. Computational cost is moderate. The principal risk is that the central temporal claim rests on one small clinical cohort.

## Portfolio B — Balanced (recommended)

1. **Kim 2018** — discovery of paired malignant residual programs and clone-aware treatment change.
2. **Yan/ARTEMIS 2026** — locked pretreatment response test plus spatial localization; never used to define programs.
3. **Wu 2021** — independent program recurrence and generic breast-cancer spatial context.
4. **Shiao 2024 (GSE246613)** — regimen-shift longitudinal falsification in TNBC.
5. **GSE205472 HR+ NAC** — negative-control subtype for generic chemotherapy/stress effects.

This portfolio maximizes complementary information per study: longitudinal clone evidence, high-powered response, spatial evidence, treatment-shift falsification, and subtype specificity. It spans institutions and platforms. Its main limitation is imperfect regimen/endpoint alignment, which is scientifically useful for falsification but complicates pooled effect sizes.

## Portfolio C — High confidence

Portfolio B plus Pal GSE161529, Zhang GSE169246, Bassez BioKey, PDX platinum time series, and GSE299631.

This improves breadth and mechanistic triangulation but increases access friction, duplicated immune-state analyses, harmonization burden, and researcher degrees of freedom. It should be used only after the core result survives Portfolio B.

## Decision

Use **Portfolio B**. Do not pool all cells or studies into a single atlas as the primary analysis. Analyze each study under a common, frozen specification and synthesize patient-level effects. Reserve Yan, Shiao, and HR+ NAC from all state definition and threshold selection.

| Property | Minimal | Balanced | High confidence |
|---|---:|---:|---:|
| Independent studies | 3 | 5 | 10 |
| Direct clinical longitudinal evidence | Low | Moderate | High |
| Spatial information | Moderate | High | High |
| Patient diversity | High | High | Very high |
| Redundancy | Low | Low | Moderate-high |
| Compute/access burden | Moderate | Moderate-high | Very high |
| External validation strength | Moderate | High | Very high if harmonization succeeds |
| Confounding risk | Moderate | Moderate | High from regimen/platform heterogeneity |
