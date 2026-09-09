# Phase 4.7 — Auditoria de candidatos longitudinais

## Candidato seleccionado para validação discreta

O GSE205472 é o candidato mais forte encontrado para a próxima validação
discreta: o manifesto GEO público fornece patient IDs, sample IDs e cinco pares
explícitos pré/pós-tratamento entre 8 doentes e 13 amostras. Não há datas/dias
biológicos, portanto `continuous_time=false` e a qualificação é
`DISCRETE_LONGITUDINAL_VERIFIED`.

O compartimento maligno é marcado como `PARTIAL`: existe scRNA de tumor, mas a
regra de identificação de células malignas ainda precisa de auditoria/fixação
independente. Por isso este candidato pode suportar validação de transições
discretas após essa auditoria, mas ainda não desbloqueia a Phase 5.

## Outros candidatos

NeoLetExe continua promissor mas controlado; SRP114962 e GSE246613 são apenas
fases ordenadas; ARTEMIS mistura material PDX e não prova biopsias repetidas do
mesmo doente; GSE176078 é cross-sectional; HTAN permanece sem proveniência
suficiente. PDX, bulk, sangue e rótulos de agregadores não foram convertidos em
dinâmica de células malignas.

Resultado: foi identificado um candidato discreto verificável, não um dataset de
tempo contínuo. A Phase 5 permanece bloqueada até validar o compartimento
maligno e a ingestão completa; não foram introduzidos modelos, simulação ou
recomendações clínicas.
