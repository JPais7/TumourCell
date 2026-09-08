# Data model

The hierarchy is `cell → sample → patient → cohort`. Every dataset declares `available_modalities` and `missing_modalities` among scRNA-seq, snRNA-seq, bulk RNA-seq, scATAC, spatial transcriptomics, IMC, Xenium, CNV, mutation, proteomics, clinical and imaging. Missing modalities are never imputed as observed. A latent record contains identity, time, continuous features, uncertainty, cell count, gene/modality coverage and QC.
