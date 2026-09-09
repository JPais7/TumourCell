# Phase 4.5.2 — Auditoria sample-level

Foi implementado um validador de registos individuais (`LongitudinalSampleRecord`), não um classificador de booleanos. O pipeline verifica identidade de amostra/biopsia, mapeamento paciente–amostra, tempos biológicos, unidades e âncoras, intervalos positivos, duplicados, tratamento e compartimento maligno.

Nesta execução não foi obtida uma tabela pública de registos sample-level verificável para SRP114962, e a metadata NeoLetExe permanece controlada no EGA. Assim, não se afirma que IDs, datas ou intervalos foram directamente inspeccionados. A evidência disponível é de estudo/publicação: SRP114962 tem fases before/mid/after NAC; NeoLetExe tem três pontos de biopsia. Ambos permanecem `ORDERED_PHASE_ONLY`.

Sequencing/run dates são retidas como metadata mas nunca convertidas em tempo biológico. Fases clínicas, ciclos de tratamento e ordem de acesso também não geram `delta_t`. O validador preserva intervalos irregulares (por exemplo 7 e 14 dias), rejeita zero/negativos e não interpola tempos ausentes.

Conclusão: `NO_CONTINUOUS_TIME_DATASET_VERIFIED`. O gate temporal falha, a Phase 5 permanece bloqueada e não foram implementados modelos, simulação, causalidade ou recomendações clínicas.
