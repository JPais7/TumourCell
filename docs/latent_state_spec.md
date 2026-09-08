# Latent-state specification

The versioned encoder is the sole new representation API. It records genes, weights, normalization, transformation, reference, missingness rules, version and a SHA-256 hash. Continuous dimensions include proliferation, IFN, antigen presentation, hypoxia, EMT, stress, stemness, inflammatory, metabolic and DNA-damage response. P8 and State_01–04 remain candidates. Tumour state is represented as `T[t] = p(z | tumour,t)` with mean, variance, quantiles, proportions, rare-state abundance, entropy and optional mixture components.

Stress/cell-cycle analyses must report unadjusted, covariate-adjusted and residualized results. Unadjusted is primary; any signal removed by adjustment is reported explicitly as potentially biological rather than silently discarded.
