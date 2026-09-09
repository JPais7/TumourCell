# TumourCell — Phase 5 dataset portfolio freeze

**Audit date:** 2026-09-10  
**Decision standard:** conservative identifiability, reproducibility and falsifiability.  
**Status:** Phase 5 can start only as a *partial* state/observation-model phase. A clone-resolved human dynamical twin is **not yet identifiable**.

## Executive conclusion

The smallest defensible portfolio is:

| Dataset | Role | Keep? | Confidence | Main reason |
|---|---|---:|---|---|
| SRP114962 / PRJNA396019 | PRIMARY HUMAN TRAINING | Yes, restricted | Medium | Human TNBC, pre/mid/post treatment and scRNA/scDNA study design; public cell-level RNA–DNA–patient linkage remains incomplete. |
| EGAS00001007242 | PRIMARY MECHANISTIC VALIDATION | Yes, controlled | High for mechanism; medium for reproducibility | Serial TNBC PDX, scRNA time series and matched single-cell CN information; strongest selection/plasticity dataset. |
| HTAPP MBC / SCP2702 | EXTERNAL REPLICATION | Yes, holdout | High | Independent human malignant-cell and TME/spatial annotations; not a treatment time series. |
| Serial spatial metastatic breast cancer dataset (2025–2026) | SPATIAL VALIDATION | Yes, guardrail | Low–Medium | Serial spatial transcriptomics with matched genomic information; only one TNBC patient in the currently identified study and preprint-level maturity. |
| GSE228154 / GSE228382 | EXPLORATORY ONLY | Defer | Medium | Public afatinib barcode/scRNA time series; lineage dynamics, not genomic clone dynamics. |
| GSE291678 / GSE291679 | EXPLORATORY ONLY | Defer | Medium–Low | Doxorubicin/ClonMapper design is useful, but a reproducible per-cell RNA–barcode table has not been reconstructed. |

Removing EGAS weakens mechanistic selection/plasticity testing. Removing SRP weakens human treatment-state training. Removing HTAPP removes independent human validation. Removing the spatial study makes the environmental/spatial component untestable. The two lineage datasets are valuable sensitivity analyses but redundant for the first core model and do not solve the genomic-clone gap.

## Evidence policy

`VERIFIED` means directly supported by a primary record, paper, or released processed object. `PARTIAL` means the study design supports the capability but the public/reproducible linkage is incomplete. `NOT_VERIFIED` means plausible but not demonstrated in the available artefacts. `NOT_AVAILABLE` means the capability is absent.

The paper-level convention that **OPcell is operative/post-treatment** is retained for the Kim curation. It must not be confused with a public ENA guarantee of patient, biopsy or continuous-time linkage; those fields remain incomplete in the archive metadata.

## Dataset audit

### 1. GSE228154 / GSE228382

**System.** TNBC preclinical barcode model/cell population, afatinib exposure, scRNA-seq and cellular barcodes. GSE228154 contains an afatinib time series (Day 0, 3, 6 and 9, with Day 0 technical replicates); GSE228382 is the associated superseries. The public GEO record exposes processed barcode/count objects. There are no human patients; the biological replicate is the experimental culture/model, not a patient.

**What is observed.** Malignant/model-cell identity is effectively known at the model level; RNA state and treatment ordering are observed; barcode lineage abundance can be compared over time. A barcode is an experimental lineage label, not a somatic genomic clone. There is no demonstrated barcode-to-genome or RNA-cell-to-DNA-cell linkage, no patient identity, no spatial context, and no patient-level biological replication.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | VERIFIED | Model is tumour-derived; no patient malignant-cell classifier is required. |
| cellular state | VERIFIED | scRNA programs and continuous scores are available. |
| state continuity | PARTIAL | Time ordering exists, but barcode/cell identity does not prove cellular continuity. |
| temporal ordering | VERIFIED | Day 0/3/6/9 treatment series. |
| treatment exposure | VERIFIED | Afatinib schedule is part of the experiment. |
| patient identity | NOT_AVAILABLE | No patient cohort. |
| clone identity | NOT_VERIFIED | Barcode is not demonstrated genomic clone identity. |
| lineage identity | VERIFIED | Cellular barcode lineage is the measured identity. |
| genomic clone identity | NOT_AVAILABLE | No matched somatic genotype. |
| RNA–DNA cell linkage | NOT_AVAILABLE | Not released as a validated paired table. |
| clonal expansion/contraction | PARTIAL | Barcode abundance changes are observable; genomic clonal selection is not. |
| state transition | PARTIAL | Population/state changes over ordered samples; individual-cell transition is not observed. |
| TME/environment | NOT_AVAILABLE | Culture model lacks patient TME. |
| spatial context | NOT_AVAILABLE | No spatial assay. |
| treatment response/resistance | PARTIAL | Afatinib tolerance/resistance phenotype is experimentally represented. |
| longitudinal pairing | VERIFIED | Experimental samples are time-ordered. |
| biological replication | PARTIAL | Replicates exist in the model, not across patients. |

**Selection versus plasticity.** Partially distinguishable: barcode enrichment supports lineage-based selection, while within-lineage transcriptional change supports plasticity. The comparison remains `selection + plasticity` at lineage level; it cannot establish genomic clone selection.

**Score (0–5).** Malignant state 3; state representation 4; temporal information 4; treatment dynamics 4; patient identity 0; clone information 2; selection identifiability 2; plasticity identifiability 3; environment/TME 0; spatial 0; human relevance 1; mechanistic value 4; external validation 1; reproducibility 4. Non-zero scores reflect public processed data and controlled time series, not human or genomic-clone evidence.

**Missing for the three equations.** No validated genomic `c`; no RNA–DNA linkage; no patient/TME `e`; no individual-cell transitions; no untreated matched control sufficient to separate adaptation from time/culture effects.

### 2. GSE291678 / GSE291679

**System.** TNBC preclinical model with ClonMapper lineage tracing, scRNA-seq and doxorubicin treatment phases. Existing repository reconstruction verifies the design and public study-level artefacts, but not a complete normalized per-cell RNA-to-barcode table.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | VERIFIED | Preclinical tumour model. |
| cellular state | VERIFIED | scRNA is available. |
| state continuity | PARTIAL | Ordered treatment phases, but no validated cell continuity. |
| temporal ordering | VERIFIED | Pre/post doxorubicin design is reported. |
| treatment exposure | VERIFIED | Doxorubicin condition is recorded. |
| patient identity | NOT_AVAILABLE | No patient cohort. |
| clone identity | NOT_VERIFIED | ClonMapper lineage is not genomic clone by default. |
| lineage identity | PARTIAL | Barcode design is reported; reconstructed per-cell linkage is incomplete. |
| genomic clone identity | NOT_AVAILABLE | No demonstrated matched somatic clone table. |
| RNA–DNA cell linkage | NOT_AVAILABLE | Not reproducibly released in the audited artefacts. |
| clonal expansion/contraction | PARTIAL | Barcode-level abundance may be compared after reconstruction. |
| state transition | PARTIAL | Population-level state changes are possible; true cell transitions are not observed. |
| TME/environment | NOT_AVAILABLE | No patient TME. |
| spatial context | NOT_AVAILABLE | No spatial assay. |
| treatment response/resistance | PARTIAL | Doxorubicin response design is present. |
| longitudinal pairing | PARTIAL | Treatment phases exist, but sample-level reconciliation is incomplete. |
| biological replication | PARTIAL | Experimental replicates, not human biological replication. |

**Selection versus plasticity.** Partially distinguishable only after a validated barcode-to-cell matrix is obtained. Even then, the result would be lineage-based clonal dynamics, not genomic evolution.

**Score.** Malignant state 3; state representation 4; temporal information 4; treatment dynamics 4; patient identity 0; clone information 2; selection identifiability 3; plasticity identifiability 3; environment/TME 0; spatial 0; human relevance 1; mechanistic value 4; external validation 1; reproducibility 2. The reproducibility score is reduced because the public reconstruction is incomplete.

**Missing.** Complete barcode–cell table, barcode QC, matched untreated controls, genomic identity, RNA–DNA linkage, and a validated sample/timepoint manifest.

### 3. SRP114962 / PRJNA396019 (Kim TNBC)

**System.** Human TNBC neoadjuvant chemotherapy study, pre/mid/post treatment design, scRNA/snRNA and scDNA components. The study reports detailed cases including eight patients and clonal extinction/persistence. The paper and supplementary tables are biologically informative; ENA/SRA metadata are not sufficient on their own to reconstruct every patient/biopsy/time relationship. The local Kim pilot contains a small processed matrix and conservative metadata, but unresolved clone identity remains unresolved.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | PARTIAL | Tumour study and malignant profiles are present, but cell-level malignant calls require the curated analysis. |
| cellular state | VERIFIED | scRNA/snRNA support program/state representation. |
| state continuity | PARTIAL | Repeated treatment phases exist; same-cell continuity is not observed. |
| temporal ordering | PARTIAL | Paper defines pre/mid/post and OPcell as post-operative; archive metadata do not fully validate biological intervals. |
| treatment exposure | PARTIAL | Neoadjuvant chemotherapy is documented at study/sample level. |
| patient identity | PARTIAL | Paper tables identify patients; public run-level reconciliation is incomplete. |
| clone identity | PARTIAL | scDNA/clonal outcomes are reported, but public cell-level RNA–clone linkage is not verified. |
| lineage identity | NOT_AVAILABLE | No experimental lineage barcode is established in the public artefacts. |
| genomic clone identity | PARTIAL | Genomic/clonal information exists at study or sample level, not a verified RNA-cell mapping. |
| RNA–DNA cell linkage | NOT_VERIFIED | This is the principal unresolved limitation. |
| clonal expansion/contraction | PARTIAL | Extinction/persistence is reported, but not as a reproducible per-RNA-cell transition table. |
| state transition | PARTIAL | Paired phase-level state change can be estimated for selected patients. |
| TME/environment | PARTIAL | RNA captures some TME context; no direct spatial ecology. |
| spatial context | NOT_AVAILABLE | No spatial transcriptomics. |
| treatment response/resistance | VERIFIED | Neoadjuvant treatment and response/resistance phenotypes are central to the study. |
| longitudinal pairing | PARTIAL | Pre/mid/post design is reported; public identifiers do not guarantee all pairings. |
| biological replication | VERIFIED | Multiple patients and samples are reported, though public reconciliation is incomplete. |

**Selection versus plasticity.** Partially distinguishable. Patient-level clonal persistence/extinction supports a selection hypothesis; paired RNA state shifts support plasticity. Without the missing RNA–DNA cell linkage, the data cannot decide whether the same genomic clone changed state or a different clone replaced it.

**Score.** Malignant state 3; state representation 4; temporal information 3; treatment dynamics 4; patient identity 2; clone information 2; selection identifiability 2; plasticity identifiability 3; environment/TME 2; spatial 0; human relevance 5; mechanistic value 3; external validation 3; reproducibility 3. The score is intentionally capped by the public linkage problem.

**Missing.** A public patient/sample/timepoint manifest linked to raw runs; validated RNA-cell ↔ DNA-cell clone mapping; exact treatment dates and dose; untreated/normal matched controls; complete malignant-cell calls; technical replicate metadata; and a reproducible small matrix with stable cell IDs.

### 4. EGAS00001007242

**System.** TNBC patient-derived xenograft, serial passaging over 2.5 years, platinum exposure, 53,641 filtered scRNA cells and matched genomic single-cell CN data from the same samples. The EGA record reports 47 controlled-access datasets. The paper explicitly analyses high-fitness clonal sweeps, weaker-fitness transcriptional plasticity, withdrawal, and hysteresis.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | VERIFIED | PDX tumour cells with genomic/CN context. |
| cellular state | VERIFIED | scRNA time series and program-level analyses. |
| state continuity | PARTIAL | Serial samples and pseudotime; not direct physical cell tracking. |
| temporal ordering | VERIFIED | Serial passaging and treatment/withdrawal design. |
| treatment exposure | VERIFIED | Platinum intervention is explicit. |
| patient identity | PARTIAL | PDX donor identity is present, but not a multi-patient human cohort. |
| clone identity | VERIFIED | Genomic/CN clone fitness is explicitly analysed. |
| lineage identity | NOT_AVAILABLE | No barcode lineage is required or established. |
| genomic clone identity | VERIFIED | Matched single-cell CN information. |
| RNA–DNA cell linkage | PARTIAL | Study reports co-measured/matched information; controlled raw access prevents independent full verification here. |
| clonal expansion/contraction | VERIFIED | Clonal sweeps and fitness differences are central results. |
| state transition | PARTIAL | State reprogramming and hysteresis are observed at serial population/cell distributions. |
| TME/environment | PARTIAL | PDX environment is controlled but not human TME. |
| spatial context | NOT_AVAILABLE | No spatial transcriptomics. |
| treatment response/resistance | VERIFIED | Platinum resistance and withdrawal are directly studied. |
| longitudinal pairing | VERIFIED | Serial PDX sampling. |
| biological replication | PARTIAL | Multiple PDX series/animals, but controlled metadata need verification. |

**Selection versus plasticity.** This is the strongest dataset for the comparison. It supports `selection + plasticity`, with genotype-associated high-fitness sweeps and CN-independent transcriptional reprogramming. It still does not establish that the same mechanism operates in primary human tumours.

**Score.** Malignant state 4; state representation 4; temporal information 5; treatment dynamics 4; patient identity 2; clone information 4; selection identifiability 4; plasticity identifiability 4; environment/TME 2; spatial 0; human relevance 1; mechanistic value 5; external validation 2; reproducibility 2. Controlled access lowers reproducibility, not the biological value.

**Missing.** Public raw access or a fully auditable derived cell-level clone/state table; human primary validation; spatial/TME context; independent untreated matched control; and explicit patient-equivalent replication.

### 5. Serial spatial metastatic breast-cancer dataset (2025–2026)

**System.** The strongest currently identified serial spatial study reports four metastatic breast-cancer patients, ten biopsies over up to 3.5 years, serial treatment timelines, CosMx 960-gene spatial transcriptomics, and patient-matched WES/CNA; only one patient is TNBC. It is currently treated as a recent/preprint-level resource until accession, peer review and all processed files are independently checked.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | PARTIAL | Tumour epithelial/CNA-supported compartments are available; exact classifier must be audited. |
| cellular state | VERIFIED | Spatial gene-expression programs and domains. |
| state continuity | PARTIAL | Serial biopsies, not same-cell tracking. |
| temporal ordering | VERIFIED | Longitudinal treatment timeline. |
| treatment exposure | VERIFIED | Treatment timelines are reported. |
| patient identity | VERIFIED | Four named longitudinal patient series. |
| clone identity | PARTIAL | Matched WES/CNA provides genomic context, not necessarily single-cell clone IDs. |
| lineage identity | NOT_AVAILABLE | No barcode lineage. |
| genomic clone identity | PARTIAL | CNA/WES support clone hypotheses; direct spatial cell clone labels require validation. |
| RNA–DNA cell linkage | PARTIAL | Patient-matched, not necessarily cell-level. |
| clonal expansion/contraction | PARTIAL | Spatial/CNA domain changes can suggest expansion; direct clonal trajectories are not guaranteed. |
| state transition | PARTIAL | Serial state distributions, not physical cell transitions. |
| TME/environment | VERIFIED | Spatial domains and neighbourhoods are directly observed. |
| spatial context | VERIFIED | Core assay capability. |
| treatment response/resistance | PARTIAL | Longitudinal resistance routes are studied; TNBC generality is limited. |
| longitudinal pairing | VERIFIED | Serial patient biopsies. |
| biological replication | PARTIAL | Four patients, one TNBC, heterogeneous disease. |

**Selection versus plasticity.** Partially distinguishable only when spatial domains and CNA are jointly analysed. It is not sufficient to identify genomic selection versus within-cell plasticity in TNBC by itself.

**Score.** Malignant state 4; state representation 4; temporal information 4; treatment dynamics 4; patient identity 2; clone information 3; selection identifiability 2; plasticity identifiability 3; environment/TME 5; spatial 5; human relevance 4; mechanistic value 3; external validation 3; reproducibility 2.

**Missing.** More TNBC patients, direct single-cell genomic linkage, full accession/processed-data verification, matched untreated controls, and adequate replication for subtype-specific inference.

### 6. HTAPP breast-cancer single-cell datasets / SCP2702

**System.** Independent human metastatic breast-cancer cohort with matched single-cell/nucleus and spatial bundles. The portal exposes counts, annotations and AnnData objects; annotations include replicate, condition, cell type, compartments, `cnv_pass_mal`, phase/QC and spatial coordinates. Existing repository audit records approximately 67 biopsies from 60 patients, with scRNA, snRNA and spatial subsets.

**Dynamic observability.**

| Capability | Status | Conservative interpretation |
|---|---|---|
| malignant identity | VERIFIED | Malignant CNV/annotation fields are available, subject to re-QC. |
| cellular state | VERIFIED | Human malignant-cell state and TME representation. |
| state continuity | NOT_AVAILABLE | No longitudinal same-patient treatment series for the intended use. |
| temporal ordering | NOT_VERIFIED | A phase field is not equivalent to validated treatment time. |
| treatment exposure | PARTIAL | Condition fields exist, but treatment-dynamic reconstruction is not established. |
| patient identity | VERIFIED | Independent human patient/sample annotations. |
| clone identity | NOT_VERIFIED | CNV/malignancy is not a longitudinal genomic clone trajectory. |
| lineage identity | NOT_AVAILABLE | No lineage tracing. |
| genomic clone identity | PARTIAL | CNV can support exploratory clone-like structure, not demonstrated clone truth. |
| RNA–DNA cell linkage | NOT_AVAILABLE | No matched single-cell DNA linkage for the core holdout. |
| clonal expansion/contraction | NOT_AVAILABLE | No validated longitudinal observation. |
| state transition | NOT_AVAILABLE | Cross-sectional validation only. |
| TME/environment | VERIFIED | Compartment and spatial data. |
| spatial context | VERIFIED | Spatial bundles and coordinates exist for subsets. |
| treatment response/resistance | PARTIAL | Useful for phenotype/state generalisation, not causal treatment dynamics. |
| longitudinal pairing | NOT_AVAILABLE | Not a validated longitudinal cohort for this model. |
| biological replication | VERIFIED | Many independent patients and technical modalities. |

**Selection versus plasticity.** Not distinguishable. HTAPP is an external observation/state/TME holdout, not a dynamic causal dataset.

**Score.** Malignant state 5; state representation 4; temporal information 1; treatment dynamics 1; patient identity 4; clone information 1; selection identifiability 1; plasticity identifiability 1; environment/TME 5; spatial 4; human relevance 5; mechanistic value 2; external validation 5; reproducibility 4.

**Missing.** Repeated treatment biopsies, exact intervention timing, matched genomic clone measurements, cell-level RNA–DNA linkage and longitudinal response labels.

## State representation

All datasets should use `RNA → continuous biological state vector s`, with programs/NMF scores as interpretable coordinates. `State_01`, `State_02`, etc. are regions of state space, not assumed discrete biological entities. Continuous scores are usable in SRP114962, EGAS, GSE228154/382, GSE291678/679 and HTAPP. Spatial data adds neighbourhood/context coordinates. RNA velocity is not treated as evidence of treatment-driven state transition unless library preparation and time scale make the assumptions defensible; no dataset receives a velocity-based claim in this audit.

## Identifiability of the three transition equations

| Equation | Best dataset(s) | Observable | Latent/non-identifiable | Verdict |
|---|---|---|---|---|
| `c_(t+1) ~ P(c_(t+1) | c_t,s_t,u_t)` | EGAS; weaker SRP; lineage datasets only as sensitivity analyses | EGAS CN clone/fitness and treatment time; SRP cohort-level persistence/extinction | Human cell-level genomic clone, migration, sampling and death processes | Partially identifiable in PDX; not identifiable in human public data. |
| `s_(t+1) ~ P(s_(t+1) | s_t,c_t,e_t,u_t)` | EGAS, SRP, GSE228/291 | RNA programs, treatment and ordered samples | Same-cell state transition, genomic clone assignment and complete environment | Partially identifiable as distributional state change. |
| `e_(t+1) ~ P(e_(t+1) | e_t,s_t,c_t,u_t)` | Spatial dataset; HTAPP for cross-sectional observation | Spatial neighbourhood/TME and serial biopsies | Environment transition, cell movement, causal niche effects | Not identifiable as a dynamic ecological equation; only partially constrained. |

The observation model `y_t ~ P_phi(y_t | c_t,s_t,e_t,m_t)` is best supported by HTAPP and the multi-modality datasets because batch, modality, depth, patient and replicate fields are available. It remains essential to model sample preparation and modality before interpreting a phase difference as biology.

## Selection versus plasticity model comparison

| Model | Evidence contribution |
|---|---|
| M1 selection only | EGAS genomic-fitness sweeps; lineage enrichment in GSE228/291; SRP persistence/extinction. |
| M2 plasticity only | EGAS withdrawal/reprogramming; SRP paired state shifts; GSE228/291 within-lineage state change. |
| M3 selection + plasticity | Best tested in EGAS; partially tested in SRP; lineage datasets support only lineage-level version. |
| M4 selection + plasticity + environment | Spatial dataset and HTAPP constrain environment observation; no public dataset fully identifies the dynamic causal environment term. |

The portfolio can compare these models, but cannot claim that M4 is identified. The correct output is a model comparison with uncertainty, not a forced mechanistic conclusion.

## Leakage and validation rules

- GSE246613/GSE176078 and any datasets already used for state discovery must not be used as external validation.
- HTAPP must remain a patient-level holdout and cannot contribute to encoder training, feature selection, hyperparameter tuning or model selection.
- Spatial data must be held out from state-definition decisions if it is used to validate environment/spatial predictions.
- Split by patient, never by cells. A cell-level random split would leak patient, sample, state and clone structure.
- Treatment phase, clone, patient and state definitions must be frozen before validation.
- All normalization, gene filtering, program discovery and malignant-cell classifiers must be fitted inside the training partition.

## Smallest Phase 5 portfolio

### Keep

1. **SRP114962** for human state/treatment training, explicitly without a claim of cell-level genomic clone transitions.
2. **EGAS00001007242** for mechanistic selection/plasticity validation.
3. **HTAPP** as an untouched external human state/TME holdout.
4. **Serial spatial dataset** only for environment/spatial validation; do not use it to train the TNBC state encoder.

### Defer

- **GSE228154/382:** keep as lineage-based sensitivity analysis for afatinib and barcode observability.
- **GSE291678/679:** keep as lineage-based sensitivity analysis for doxorubicin, pending a reproducible per-cell barcode matrix.

They are not included in the core because both are preclinical lineage datasets and are redundant with EGAS for the first mechanistic model. They can later test whether conclusions depend on treatment or on lineage-label technology.

## Architecture support matrix

| Dataset | State | Clone | Treatment | Time | Selection | Plasticity | Environment | Spatial | Human |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| GSE228154/382 | 4 | 2 | 4 | 4 | 2 | 3 | 0 | 0 | 1 |
| GSE291678/679 | 4 | 2 | 4 | 4 | 3 | 3 | 0 | 0 | 1 |
| SRP114962 | 4 | 2 | 4 | 3 | 2 | 3 | 2 | 0 | 5 |
| EGAS00001007242 | 4 | 4 | 4 | 5 | 4 | 4 | 2 | 0 | 1 |
| Serial spatial 2025–2026 | 4 | 3 | 4 | 4 | 2 | 3 | 5 | 5 | 4 |
| HTAPP/SCP2702 | 4 | 1 | 1 | 1 | 1 | 1 | 5 | 4 | 5 |

| Architecture component | Dataset constraints |
|---|---|
| A. State encoder `RNA → s` | SRP training; EGAS mechanistic state; HTAPP external holdout. |
| B. Evolution `c_(t+1) | c_t,s_t,u_t` | EGAS primarily; SRP only at partially observed cohort level. |
| C. Plasticity `s_(t+1) | s_t,c_t,e_t,u_t` | EGAS and SRP; GSE228/291 sensitivity analyses. |
| D. Ecological `e_(t+1) | e_t,c_t,s_t,u_t` | Spatial dataset constrains observations; dynamic causal term remains unresolved. |
| E. Observation `y_t | c_t,s_t,e_t,m_t` | HTAPP and multi-modality metadata support batch/patient/modality effects. |

## Gap analysis

**Can the full Phase 5 start? No.** A preliminary Phase 5 state-space prototype can start, but the stated clone-resolved human twin cannot be claimed as identifiable.

Remaining biological bottlenecks:

1. No public human TNBC dataset with reproducible cell-level RNA ↔ genomic-clone linkage across repeated treatment timepoints.
2. No direct observation of individual-cell state transitions; serial samples provide distributions, not trajectories.
3. The environment transition equation is under-observed; spatial studies are sparse and subtype-limited.
4. Treatment dose, exact dates, sampling intervals and technical replicate structure are incomplete for SRP114962.
5. PDX mechanistic evidence does not establish human causal equivalence.
6. Selection and plasticity remain confounded in human data by sampling, death, migration and clone misassignment.

## Stopping criterion and falsifiable prediction

The portfolio is sufficient to begin a restricted Phase 5 when it produces one prediction that can be wrong:

> After patient-level and library-depth adjustment, the pre-treatment malignant state distribution in SRP114962 predicts a treatment-associated shift in the mid/post state distribution; in EGAS, the direction and persistence of that shift differs between high-CN-fitness and weak-fitness clones.

**Falsification:** the state score does not predict the later distribution in held-out patients, or EGAS shows no clone-fitness-dependent difference after treatment/withdrawal, or the direction reverses under a pre-registered analysis.

**Test dataset:** SRP114962, with patient-level cross-validation.  
**Mechanistic validation:** EGAS00001007242.  
**Independent state/TME validation:** HTAPP.  
**Spatial guardrail:** serial spatial dataset.

## Most important next experiment

A prospective TNBC cohort with pre/mid/post biopsies, joint scRNA + scDNA or validated single-cell multi-omics from the same aliquot, explicit treatment dates/doses, patient/sample/barcode reconciliation, matched spatial transcriptomics, and technical replicates. This single experiment would reduce the dominant uncertainty: whether a later RNA state belongs to the same genomic clone, a different clone, or a changed environment.

## Sources

- [GSE228154 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228154) and [GSE228382 superseries](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE228382)
- [GSE291678 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291678) and [GSE291679 GEO record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE291679)
- [Kim study PRJNA396019](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA396019)
- [EGAS00001007242 EGA record](https://ega-archive.org/studies/EGAS00001007242)
- [HTAPP MBC / SCP2702](https://singlecell.broadinstitute.org/single_cell/study/SCP2702/htapp-mbc)
- [Serial spatial metastatic breast-cancer study](https://pmc.ncbi.nlm.nih.gov/articles/PMC13105127/)
- Local supporting audits: `docs/phase4_5_5_authoritative_srp114962_resolution.md`, `docs/phase4_9_5_cell_to_genome_linkage.md`, `docs/phase4_9_12_clonmapper_reconstruction.md`, `docs/phase3c_dataset_audit.md`.

