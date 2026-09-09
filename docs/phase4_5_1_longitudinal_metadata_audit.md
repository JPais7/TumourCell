# Phase 4.5.1 — Auditoria de metadata longitudinal

## Conclusão

Não foi verificado nenhum dataset candidato com `days_from_treatment` ou datas de biopsia públicas, ao nível paciente/amostra, que permitam calcular intervalos biológicos. “Before/mid/after” e “three timepoints” demonstram ordem clínica, não intervalos iguais. Datas de sequenciação/run não são usadas como tempo biológico.

| Dataset | Repetição | Metadata temporal | Tempo biológico verificado | Tratamento | Maligno | Clone | Classificação |
|---|---|---|---|---|---|---|---|
| SRP114962 | Sim; 20 doentes, seis explicitamente no snRNA | before/mid/after NAC | Não | Sim | Sim | scDNA/scRNA | ORDERED_PHASE_ONLY |
| NeoLetExe / EGAS50000001362 | Sim; 25 doentes, três biopsias | três pontos durante inibidor de aromatase | Não; metadata controlada | Sim | Sim | scTCR/scBCR | ORDERED_PHASE_ONLY |

SRP114962 é suportado pelo registo SRA/PRJNA396019 e pela descrição do estudo: amostras longitudinais durante quimioterapia, mas sem uma tabela pública verificada de dias por biopsia. A página EGA e a documentação OUH confirmam 25 doentes e três pontos de biopsia no NeoLetExe, mas a metadata detalhada exige acesso controlado. A classificação não é promovida por inferência.

`ACCESS_CONTROLLED_METADATA_REQUIRED` distingue “não verificado” de “não existente”. O próximo passo é obter attributes BioSample/ENA e metadata EGA autorizada. Mesmo que sejam encontradas datas, não se inicia forecasting nesta fase: primeiro é necessário validar identidade, tratamento, malignidade, duplicados e intervalos desiguais.

O `TEMPORAL_DATASET_GATE` permanece **NOT_PASSED** e a Phase 5 permanece bloqueada. Não foi criada dinâmica contínua, simulação, inferência causal ou recomendação clínica.
