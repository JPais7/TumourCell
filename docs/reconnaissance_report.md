# Scientific reconnaissance report

Date: 2026-09-06

## Executive decision

The broad central question is **partly answerable** with public data. Public cohorts can establish reproducible cellular programs, compare pretreatment and post-treatment samples, and test external prediction. They cannot, in general, recover patient-specific continuous dynamics or prove cell-state transitions because most clinical studies contain only two or three destructive biopsies and no lineage tracing.

TNBC is the strongest initial context because it has the best combination of biological importance, neoadjuvant response endpoints, paired single-cell studies, spatial cohorts, and independent validation opportunities. This is a data-driven choice, not a claim that TNBC biology is intrinsically more modelable than ER-positive disease.

However, a broad TNBC state/response paper would now be incremental. Yan et al. profiled 427,857 cells from 101 treatment-naive TNBC patients, with spatial data in 44, and linked malignant metaprograms and multicellular ecotypes to neoadjuvant response ([Nature 2026](https://www.nature.com/articles/s41586-026-10469-9)). Kim et al. already reported that resistant genotypes pre-existed while transcriptional phenotypes were acquired under chemotherapy in longitudinal TNBC ([Cell 2018](https://pubmed.ncbi.nlm.nih.gov/29681456/)). The defensible gap is therefore a cross-study, patient-level falsification of the **relative contribution and reproducibility of selection versus transcriptional plasticity**.

## Answers to the checkpoint questions

1. **Answerable with public data?** Yes for reproducible states, treatment-associated changes, and external prediction. Only partially for transitions; no for causal transition rates or a clinical digital twin.
2. **Strongest subtype?** TNBC. ER+/HER2− is the strongest alternative because FELINE provides serial endocrine/CDK4/6 biopsies, but public accessibility and independent single-cell validation are weaker and treatment mechanisms differ.
3. **Most relevant studies?** Fifteen prioritized studies are recorded in `dataset_utility_matrix.csv`: Yan/ARTEMIS, Kim, Wu, Pal, Karaayvaz, Bassez, Zhang, Shiao, Wang, Hammerl, Bassiouni, HR+ NAC GSE205472, FELINE GSE114727, TNBC PDX GSE276609/associated series, and the BRCA1-deficient residual-disease atlas GSE299631.
4. **Complementary studies?** Yan supplies high-powered pretreatment response/spatial context; Kim supplies paired human chemotherapy and clone-aware RNA/DNA; Wu and Pal supply cross-subtype state recurrence; Shiao and GSE169246 supply longitudinal immunotherapy contexts; PDX time series supplies dense perturbation timing.
5. **Redundant studies?** Wu and Pal overlap as untreated atlases but differ in sampling and scale. Bassez and Zhang both address checkpoint therapy but are not interchangeable: treatment, assay, cohort, and access differ. Downstream reanalyses of GSE176078 are not independent validation.
6. **Minimum sufficient portfolio?** Kim + Yan + Wu + one locked independent treatment cohort (Shiao for immunotherapy generality or GSE205472 for subtype negative control). See `minimum_sufficient_portfolio.md`.
7. **Reserve for external validation?** Yan/ARTEMIS must be held out if Kim defines residual programs; alternatively reverse these roles, but never tune on both. Shiao GSE246613 is a stringent regimen-shift validation. Pal GSE161529 is a state-recurrence validation, not response validation.
8. **Useful longitudinal information?** Kim (pre/mid/surgery), Bassez (pre/on anti-PD1), Zhang GSE169246 (pre/post paclitaxel ± atezolizumab), Shiao GSE246613 (pre/after pembrolizumab/after pembrolizumab+radiation), FELINE GSE114727 (serial endocrine ± ribociclib), GSE205472 (paired HR+ NAC), and GSE299631 (residual disease; verify human pairing).
9. **Spatial information?** Yan/ARTEMIS (Xenium, Visium, Visium HD), Wu (Visium), Shiao (spatial proteomics), Wang (multiplex spatial profiling), Hammerl (spatial immunophenotyping), Bassiouni (spatial transcriptomics), GSE299631 (ST/IMC).
10. **Strongest confounders?** Patient, treatment regimen, biopsy site, primary versus metastasis, response definition, biopsy timing, dissociation bias, tumour purity, inferred-malignant-cell errors, cell cycle, hypoxia/stress, sample viability, platform, chemistry, institution, and compositionality.
11. **Competing explanations?** Selection of a pre-existing clone; reversible stress; proliferation arrest; hypoxia; immune/stromal admixture; site-specific ecology; technical loss of fragile populations; patient-specific CNA programs; and sampling of different tumour regions.
12. **Genuinely novel?** Potentially: a preregistered, clone-aware, cross-study estimate of how often residual transcriptional programs are detectable before therapy and whether that estimate generalizes across patients and platforms. Novelty is conditional on a final systematic search and metadata audit.
13. **Already known?** TNBC is heterogeneous; malignant programs and ecotypes associate with response; resistant genotypes can pre-exist; treatment can reprogram transcriptomes; immune/spatial niches associate with immunotherapy response.
14. **What falsifies the hypothesis?** Locked residual programs fail to enrich post-treatment in an independent cohort; apparent emergence vanishes after controlling for CNA clone, stress, cell cycle, purity, or site; pretreatment detectability estimates are no better than downsampling controls; results reverse under reasonable state definitions.
15. **Kill criteria?** Listed in `kill_criteria.md`.
16. **Smallest credible v0.1?** A patient-level, clone-aware state-program reproducibility and transition-accounting analysis across one discovery and one external validation study, with a negative-control subtype or regimen.
17. **Do not attempt yet:** neural ODEs, patient-specific forecasts, causal ligand–receptor claims, clinical treatment recommendations, multi-omic end-to-end integration, or virtual drug ranking.
18. **Strongest publication question:** Are residual TNBC programs predominantly selected from detectable pretreatment states or induced within persistent clones, and is the balance reproducible across independent clinical cohorts?
19. **Strongest later virtual perturbation question:** Which program regulators are predicted—consistently across simple transition models and perturbational reference data—to reduce entry into a validated residual program without merely suppressing proliferation? This is postponed until Level 2 succeeds.
20. **Future experiment:** Barcode a diverse TNBC organoid/PDX panel, sample before/during/after chemotherapy, jointly measure lineage and transcriptome, and perturb a locked regulator. Test whether the same clone changes state and whether blocking the regulator lowers residual-state occupancy without nonspecific cytotoxicity.

## Recommended scientific specification

**Central question.** Across independent TNBC cohorts, do malignant-cell programs enriched in residual disease arise mainly by expansion of detectable pretreatment programs, or by treatment-associated transcriptional reprogramming within persistent CNA-defined clones?

**Primary hypothesis.** Both mechanisms occur, but a reproducible stress-adjusted residual program appears within persistent genomic clones and cannot be explained solely by pretreatment abundance, cell cycle, hypoxia, or sampling depth.

**Discovery strategy.** Use Kim et al. paired samples to define within-patient malignant programs with NMF/cNMF-like program discovery. Infer broad CNA clones cautiously; quantify program usage rather than hard clusters. Model treatment, time, and clone with patient-blocked statistics. Use downsampling to estimate whether “absent before treatment” is distinguishable from “undetected.”

**Validation strategy.** Freeze genes, signs, scoring, QC, and thresholds. Test post-treatment enrichment and pretreatment detectability in an independent cohort at the patient level. Use Yan only for pretreatment response association and spatial/ecotype context, not as evidence of temporal emergence. Use leave-one-patient-out internal checks only as secondary evidence.

**Falsification strategy.** Predefine nulls: label permutation within study constraints, matched cell-depth downsampling, generic stress/cell-cycle scores, randomized gene sets matched for expression, and alternative malignant-cell/CNA calls. Require directionally consistent patient-level effects in an independent study.

**v0.1 architecture.** Immutable raw/metadata manifests → study-specific QC → malignant-cell and CNA evidence → interpretable program scores → patient/sample pseudobulk → discrete transition accounting → locked external test → sensitivity/falsification report. No continuous-time simulator.

**Phase 1 plan.** (1) Verify accessions, consent/access, patient/sample maps, and raw layers. (2) Preregister estimand, split, state definition, and falsifiers. (3) Download only Kim plus the selected independent validation cohort. (4) Reproduce published cohort summaries. (5) Run program robustness and confounder tests. (6) Stop at the decision gate and report positive and negative findings.

## Subtype comparison

| Criterion | TNBC | ER+/HER2− | HER2+ |
|---|---|---|---|
| Public single-cell breadth | Strong | Moderate | Limited–moderate |
| Serial treatment tissue | Several small/medium cohorts | FELINE and HR+ NAC, but fewer independent public cohorts | Sparse |
| Spatial response data | Strong and rapidly expanding | Limited | Limited |
| Biological endpoint | pCR/residual disease is clear | Ki67/PEPI/resistance definitions vary | Response confounded by targeted regimens |
| Novelty room | Narrow after 2026 ecotype study | Potentially wider | Potentially wider but data-limited |
| Initial decision | **Selected** | Strong alternative/negative control | Not selected |

## Gate assessment

| Gate | Status | Reason |
|---|---|---|
| Data | Conditional pass | Multiple public processed datasets; raw/metadata access must be audited. |
| Novelty | Conditional pass | Broad state-response analysis fails; clone-aware cross-study falsification may pass. |
| Identifiability | Partial pass | Discrete observed changes are identifiable; continuous transition rates are not. |
| Validation | Pass in principle | Independent cohorts/platforms exist, though exact endpoint alignment is imperfect. |
| Falsification | Pass | Strong technical and biological nulls are available. |
| Scientific value | Pass | Distinguishing selection from plasticity would change how residual-state biomarkers and perturbations are interpreted. |

## Strongest rejected direction

**Rejected:** build a pretreatment TNBC ecotype predictor of chemotherapy response. It has excellent data and immediate clinical relevance, but Yan et al. 2026 already defines malignant metaprograms, TME states, ecotypes, spatial organization, and a response panel in a much larger cohort than a v0.1 could assemble. Reanalysis would be replication unless it asks a materially different, externally validated question.

## Evidence boundaries

- “Transition” means a difference between observed sample-time distributions, not tracked movement of the same cell.
- CNA inference from RNA is supporting evidence, not ground-truth lineage.
- Patient is the primary independent unit; cells improve measurement precision, not biological sample size.
- The report is reconnaissance, not a completed systematic review. Every `verify` field must be resolved against publication supplements and repository metadata before analysis.
