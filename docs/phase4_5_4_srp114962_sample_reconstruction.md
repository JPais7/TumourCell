# Phase 4.5.4 — Reconstrução sample-to-patient do SRP114962

## Hierarquia e fontes

Foi consultada directamente a API pública ENA para `SRP114962`. O resultado
contém 6.633 runs, 6.633 experiments, 6.633 `sample_accession` e 6.633
`biosample_accession` (o accession SAMN foi usado como referência técnica, não
como prova de paciente). Foram inspeccionados registos representativos
`SAMN09909196`, `SAMN09909195` e `SAMN09909424`, com aliases
`KTN206OPcell108`, `KTN206OPcell107` e `KTN206OPcell336`.

Runs e experiments não foram confundidos com doentes ou biopsias. Os aliases
`KTN...cell...` não foram agrupados em pacientes sem campo autoritativo. A
coluna ENA `collection_date` está vazia nos 6.633 registos. Não foi encontrada
uma tabela suplementar pública de mapping sample→patient→biopsy→timepoint;
`srp114962_supplementary_mapping.csv` contém apenas o cabeçalho.

## Resultado

Não há patient IDs, biopsy IDs, timepoints ou tempos biológicos directamente
verificados. Portanto: `verified_patient_count=0`,
`patients_with_repeated_samples=0`, `direct_biological_time_records=0` e
`valid_intervals=0`. O dataset permanece `ORDERED_PHASE_ONLY`, sustentado apenas
pela descrição longitudinal da publicação, não por tempo contínuo observado.

NeoLetExe permanece `ACCESS_CONTROLLED_METADATA_REQUIRED`; não foi inferido
qualquer valor controlado. A Phase 5 continua bloqueada.
