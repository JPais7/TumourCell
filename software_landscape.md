# Software landscape

| Tool/approach | Appropriate use | Does not establish | Decision |
|---|---|---|---|
| [Scanpy](https://scanpy.readthedocs.io/en/stable/) | Reproducible preprocessing, QC, visualization, neighborhood graphs and basic trajectory tools at atlas scale | Biological independence or causal dynamics | Use as core data layer |
| [scvi-tools/scVI](https://docs.scvi-tools.org/en/stable/user_guide/models/scvi.html) | Count-aware latent representation, batch covariates, reference mapping; scales beyond one million cells | Interpretable states by itself; correct removal of confounded biology | Use selectively; benchmark against non-integrated program scoring |
| cNMF/NMF programs | Interpretable recurrent gene programs that allow mixed state usage | Lineage or direction | Preferred state representation for v0.1 |
| [scVelo](https://scvelo.readthedocs.io/en/latest/about.html) | Exploratory directionality when spliced/unspliced layers and kinetic assumptions are credible | Clinical longitudinal transitions or causality | Do not make primary; use only as sensitivity analysis |
| [CellRank](https://cellrank.readthedocs.io/en/stable/notebooks/tutorials/general/100_getting_started.html) | Fate probabilities from velocity, time, pseudotime, lineage or multiview kernels | Ground-truth fate without informative kernels | Postpone; eligible only if input evidence passes diagnostics |
| Milo/scCODA-like methods | Differential abundance with neighborhoods or compositional modeling | Within-cell reprogramming | Use patient-aware formulation where compatible |
| Squidpy/SpatialData | Spatial neighborhoods, graphs and image-linked data management | Temporal constraints from cross-sectional sections | Use for locked spatial localization |
| HMM/Markov state models | Discrete transitions when repeated observations identify transition probabilities | Patient-specific rates from two destructive biopsies | Not justified in v0.1 |
| Mechanistic ODE/agent-based models | Explicit assumptions and future perturbation hypotheses | Calibration without repeated quantitative state data | Postpone |
| [TumorTwin](https://pubmed.ncbi.nlm.nih.gov/42116079/) | Modular patient-specific oncology model initialization/update | Evidence that this biological problem is identifiable | Monitor/adapt later; do not adopt now |
| MRI-based breast twin | Tumour-scale growth/response calibration ([npj Digital Medicine 2025](https://www.nature.com/articles/s41746-025-01579-1)) | Cell-state transition dynamics | Complementary future scale, not v0.1 |

Build new code only for study manifests, locked program transfer, depth-aware nondetection, discrete selection/plasticity accounting, and falsification reports. Reuse established infrastructure elsewhere.
