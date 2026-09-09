# Phase 4.5.3 — Aquisição de metadata longitudinal real

## Aquisição directa

Foi consultada a API pública ENA para `SRP114962`. Foram examinadas 6.633
linhas `read_run`, com 6.633 `sample_accession` únicos. A coluna
`collection_date` estava vazia em todas as linhas. Foram obtidos identificadores
de amostra e aliases técnicos, mas não patient IDs, biopsy dates, treatment
anchors ou elapsed biological time.

Para NeoLetExe (`EGAS50000001362`), a descrição EGA/OUH confirma três pontos de
biopsia e 25 doentes, mas não foi obtida uma tabela pública sample-level; os
campos essenciais permanecem `ACCESS_CONTROLLED_METADATA_REQUIRED`.

## Resultado

| Dataset | Registos examinados | Tempo biológico directo | Intervalos válidos | Classificação |
|---|---:|---:|---:|---|
| SRP114962 | 6.633 | 0 | 0 | ORDERED_PHASE_ONLY |
| NeoLetExe | 0 públicos | 0 | 0 | ACCESS_CONTROLLED_METADATA_REQUIRED / ORDERED_PHASE_ONLY |

Nenhum intervalo foi calculado. Não foram usados sequencing/run dates, ordem de
acessão, ciclos de quimioterapia ou fases clínicas como `delta_t`. O ficheiro
`longitudinal_patient_intervals.csv` contém apenas o cabeçalho, porque não há
intervalos directamente observados elegíveis.

Conclusão: `NO_CONTINUOUS_TIME_DATASET_VERIFIED`; a Phase 5 continua bloqueada.
