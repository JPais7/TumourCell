# Kim 2018 supplement extraction

The first structured extraction from `mmc1.pdf` records the 20-patient clinical
cohort and the eight patients with single-cell DNA and RNA profiling. The study
labels P1, P2, P6 and P9 as clonal
extinction cases, and P11, P12, P14 and P15 as clonal persistence cases.

The three ordered treatment phases are pre-treatment (0 cycles), mid-treatment
(2 cycles), and post-treatment/surgery (6 cycles). The source supplement contains
Tables S1-S6: clinical cohort data, exome coverage, indels, purity/cellularity,
post-treatment mutations, and targeted amplicon validation.

Tables S1 and S4 have now been transcribed visually into
`data/manifests/kim_table_s1_clinical.csv` and
`data/manifests/kim_table_s4_purity.csv`. The transcription preserves `NA`,
`Not done`, and empty assay fields rather than converting them to zero.

The ENA `read_run` metadata for `SRP114962` was queried without downloading
sequence files. Matching `sample_alias` or `sample_accession` against the 20
clinical IDs in Table S1 recovered all 20 patients: 6,626 of 6,633 runs. The
eight-patient single-cell subset is represented by 6,586 runs and an
ENA-reported FASTQ byte total of approximately 118.1 GB. This is a metadata
estimate, not a download requirement. The selected subset contains RNA-seq,
WGS and WXS records; the alias convention includes patient IDs and cell
identifiers, but treatment timepoint must still be reconciled against the
paper's sample tables and cannot be inferred solely from the alias.

The complete crosswalk is in `data/manifests/kim_ena_alias_matches.csv`; the
eight-patient subset is in `data/manifests/kim_8patient_ena_metadata.csv`.
Cell-level clone identifiers and RNA-to-DNA linkage remain unresolved.

Alias-level timepoint audit: the eight-patient RNA-seq subset contains
`0cell` aliases in all eight patients, `2cell` aliases in KTN132, KTN302 and
KTN615, and `OPcell` aliases in the remaining records. No aliases explicitly
state `pre`, `mid`, or `post`. The strings `0cell` and `2cell` are compatible
with the paper's 0-cycle and 2-cycle phases, but are not promoted to formal
timepoint labels without an independent key. The post-treatment phase is not
identified by ENA aliases. The full audit is in
`data/manifests/kim_8patient_timepoint_reconciliation.csv`.

For downstream selection, a conservative labelled manifest was generated at
`data/manifests/kim_8patient_conservative_timepoints.csv`: 2,657 runs are
`pre_candidate`, 1,151 are `mid_candidate`, and 1,920 remain `unknown`.
These are candidate labels with explicit confidence/basis fields, not
asserted biological timepoints.
