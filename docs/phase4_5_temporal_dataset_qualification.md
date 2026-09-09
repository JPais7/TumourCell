# Phase 4.5 — Qualificação de datasets temporais

Esta fase avalia observabilidade; não constrói um modelo de dinâmica. `CONTINUOUS_TIME_OBSERVABLE` exige repetição no mesmo doente, identificadores de doente/amostra, tempo biológico decorrido verificável, tratamento e compartimento maligno. Uma data de sequenciação não é uma data biológica.

| Dataset | Repetição | Tempo | Tratamento | Maligno | Clone | Classificação |
|---|---:|---:|---:|---:|---:|---|
| GSE246613 | sim | fases | sim | sim | desconhecido | ORDERED_PHASE_ONLY |
| SRP114962 | sim | fases, intervalos por verificar | sim | sim | sim | ORDERED_PHASE_ONLY |
| NeoLetExe | sim | 3 fases, dias por verificar | sim | sim | sim | ORDERED_PHASE_ONLY |
| ARTEMIS | não no scRNA público | desconhecido | sim | sim | desconhecido | NOT_LONGITUDINAL |
| GSE176078 | não | não | desconhecido | sim | desconhecido | NOT_LONGITUDINAL |
| HTAN metastatic | desconhecido | desconhecido | desconhecido | desconhecido | desconhecido | UNKNOWN |

GSE246613 permanece `ORDERED_PHASE_ONLY`: Base, PD1 e RTPD1 fornecem ordem, não `delta_t`. SRP114962 e NeoLetExe são os candidatos mais promissores porque têm amostras repetidas e informação clonal/terapêutica, mas não foram promovidos sem verificar metadados de datas/intervalos. O acesso EGA controlado também não é tratado como disponibilidade imediata.

O gate `TEMPORAL_DATASET_GATE` falha actualmente. Não existe candidato verificado que cumpra simultaneamente tempo decorrido, repetição, tratamento e compartimento maligno; portanto a Phase 5 continua bloqueada.

Patient generalization (doentes separados), temporal generalization (tempo futuro) e joint patient+temporal generalization são propriedades distintas. Esta fase apenas qualifica se a terceira poderá ser tentada no futuro.

O ranking é qualitativo e transparente: `BEST_CANDIDATE` requer tempo contínuo verificado e malignidade; `PROMISING` corresponde a fases ordenadas; `PARTIAL` a longitudinalidade insuficiente; `INSUFFICIENT` aos restantes. Não há score arbitrário, não há dados sintéticos e continuam desactivados simulação, causalidade e recomendações clínicas.
