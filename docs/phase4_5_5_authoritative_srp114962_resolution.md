# Phase 4.5.5 — Resolução autoritativa do SRP114962

## Fontes e hierarquia

Foram consultados directamente a API ENA para `SRP114962` e os registos NCBI
BioSample (`SAMN...`) correspondentes aos 6.633 acessions. Todos os 6.633
BioSamples foram obtidos em lotes E-utilities e os atributos foram contados.
Também foram consultados SRA/PRJNA396019 e a descrição publicada do estudo.

Os BioSamples fornecem atributos directos `treatment` e `sample_type`; não
fornecem nos registos inspeccionados patient ID, biopsy ID, timepoint ou data
biológica. Foram observados 829 registos `sample_type=tumor`. Os valores de
treatment foram: pre 2.682, mid 1.151, post 1.920, preTX 371, midTX 275,
postTX 183, blood 21, operative 19 e 2cycleschemo 11.

`KTN...cell...`, `SRR...`, `SRX...` e `SAMN...` permanecem identificadores
técnicos. Não foram agrupados em pacientes por prefixo ou numeração. Não foi
encontrada uma tabela suplementar pública de mapping; o artefacto correspondente
contém apenas o cabeçalho. Não há datas de recolha nem âncoras de tratamento
ligadas a uma biopsia individual.

## Decisão

`verified_patient_count=0`, `patients_with_repeated_samples=0`,
`direct_biological_time_records=0`, `valid_intervals=0`. O tratamento é
observável como categoria de amostra, mas não como relógio biológico.

Classificação final: **ORDERED_PHASE_ONLY**. Não é permitido inferir paciente,
biopsia ou tempo contínuo a partir dos aliases. A Phase 5 permanece bloqueada.
