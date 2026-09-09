# Phase 4.4.1 — Semântica temporal real e baselines corrigidos

## Resultado

O benchmark corrigido do GSE246613 usa a representação congelada e preserva 34
doentes, 99 observações e 65 transições adjacentes. A fonte disponibiliza apenas
fases clínicas ordenadas (`Base → PD1 → RTPD1`), não datas nem dias decorridos.

Por isso, `time_value` é explicitamente ordinal (`ordered_clinical_phase`) e
`delta_t` é `NOT_AVAILABLE`. Não são assumidos intervalos unitários e não é
ajustada dinâmica contínua. O resultado é descritivo e não causal.

## Baseline

O baseline de persistência prevê delta zero. A sua MAE é **0.9272538899566909**
e RMSE **1.3508372243554625** na representação congelada. A acurácia
direccional é `NOT_APPLICABLE`, porque o baseline não prevê aumento ou redução.
Cobertura de intervalos 50/80/95 é `NOT_ESTIMABLE` sem um modelo de incerteza
válido.

O benchmark numérico anterior foi **superseded**: usava `delta_t = 1` para
categorias clínicas e é formalmente **INVALIDATED_BY_TEMPORAL_SEMANTICS_CORRECTION**.

O benchmark actual é um `PATIENT_HELD_OUT_ORDERED_PHASE_BENCHMARK`, não uma
validação de dinâmica temporal contínua. As etiquetas de resposta são mantidas
apenas como metadados e não entram em estados, preditores, splits ou selecção.

## Limitação e próximo requisito

Para estimar dinâmica contínua e forecasting prospectivo são necessários dias
desde a biópsia/tratamento (ou timestamps harmonizados) por amostra. Não se
inicia a Fase 5 com os dados atuais. Também serão necessários observações
repetidas suficientes, metadados de tratamento e, idealmente, resolução clonal e
espacial/contextual.
