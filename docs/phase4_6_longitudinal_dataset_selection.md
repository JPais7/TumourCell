# Phase 4.6 — Selecção de dataset longitudinal real

Foram investigados seis candidatos: NeoLetExe, SRP114962, ARTEMIS, GSE176078,
HTAN metastatic breast cancer e GSE246613. A evidência foi separada em
descrição de publicação, metadata de repositório e metadata sample-level. Um
rótulo “longitudinal” de agregador nunca foi tratado como prova.

## Decisão

Nenhum candidato cumpre actualmente simultaneamente patient linkage autoritativo,
repetição de amostras biológicas, compartimento maligno verificável e tempo
biológico/discreto ao nível da amostra. Portanto o resultado é
**NO_PUBLIC_LONGITUDINAL_DATASET_VERIFIED**.

NeoLetExe é o candidato mais promissor, mas a metadata EGA sample-level está
controlada. SRP114962 foi auditado directamente: os BioSamples têm categorias de
tratamento, mas não patient IDs, biopsias ou tempos. ARTEMIS inclui material PDX
longitudinal, que não é automaticamente uma série de biopsias do mesmo doente.
GSE176078 é cross-sectional; HTAN permanece com proveniência insuficiente;
GSE246613 é `ORDERED_PHASE_ONLY`.

O ranking é qualitativo e fail-closed: campos críticos ausentes não podem ser
compensados por grande número de células, resposta clínica ou disponibilidade de
PDX. A Phase 5 continua bloqueada; não foram iniciados modelos, simulação ou
recomendações clínicas.
