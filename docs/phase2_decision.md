# Decisão da Fase 2

Data: 2026-09-06  
Decisão global: **NÃO AVANÇAR para modelação dinâmica, virtual perturbation ou Digital Twin**

## Tabela de decisão

| Questão | Resultado | Evidência | Confiança |
| --- | --- | --- | --- |
| P8 replica em ARTEMIS? | **INCONCLUSIVO** | +0,187 SD; IC95% −0,242 a +0,607; critérios congelados não cumpridos | Moderada para o veredicto; baixa para o efeito real |
| Direção preservada? | **Sim na análise primária; não em todas as sensibilidades** | 100% LOO positivo na primária; −0,010 no threshold ≥100 e −0,002 ajustado | Moderada |
| Magnitude compatível? | **Amplamente compatível, mas imprecisa** | ARTEMIS +0,187 versus discovery +0,343 SD; IC amplos e endpoints apenas parcialmente comparáveis | Baixa |
| Resultado depende de poucos doentes? | **Não na análise primária; depende da especificação de QC** | LOO +0,143 a +0,269; P81 é o mais influente, sem inversão; threshold ≥100 elimina o efeito | Moderada |
| CNV é identificável? | **Parcialmente, apenas CNV-like de RNA (níveis 1–3)** | CopyKAT/consenso publicado; sem validação de DNA nos objetos analisados | Moderada |
| Clonalidade é identificável? | **NÃO IDENTIFICÁVEL** | Clusters RNA-CNV não equivalem a clones biologicamente validados | Alta |
| Seleção vs plasticidade é identificável? | **NÃO IDENTIFICÁVEL** | Sem rastreio clonal baseline→residual; ARTEMIS molecular é pré-tratamento | Alta |

## OBSERVADO

- P8 foi aplicado sem retuning e passou a cobertura genética pré-especificada.
- O efeito ARTEMIS primário aponta para P8 mais alto no baseline dos doentes pCR, mas o IC95% inclui zero.
- Nenhum doente isolado inverte a análise primária.
- O efeito enfraquece até aproximadamente zero no threshold de 100 células e no ajuste técnico conjunto.
- A inferência CNV disponível é baseada em expressão; genes HLA, centrais em P8, são um risco documentado de CNV focal espúrio.

## INFERIDO

- Existe um sinal independente fraco, direcionalmente compatível com a hipótese biológica, mas insuficiente para estabelecer generalização robusta.
- A estabilidade leave-one-out sugere que o sinal primário não é um artefacto de um único doente; a instabilidade de QC/técnica limita a sua transportabilidade.
- A classificação aneuploide é defensável para enriquecer células malignas, não para atribuir identidade clonal.

## SUPORTADO

- Aplicabilidade técnica da definição P8 congelada a ARTEMIS.
- Comparabilidade apenas parcial entre o endpoint baseline pCR/RD e o fenómeno longitudinal R2/NR.
- Necessidade de manter P8 como hipótese, não como biomarcador validado.

## INCERTO

- Se o efeito verdadeiro de P8 em ARTEMIS é positivo e biologicamente relevante.
- Se a atenuação ajustada representa confundimento técnico, perda de precisão/colinearidade ou ausência de associação independente.
- Se P8 generaliza a outros regimes e coortes com endpoints mais homólogos.

## NÃO IDENTIFICÁVEL

- Clones biologicamente defensáveis a partir dos objetos scRNA-seq atuais.
- Correspondência do mesmo clone entre baseline e residual.
- Seleção de um estado P8-like versus aquisição plástica de P8 após tratamento.
- Causalidade ou mecanismo terapêutico.

## Aplicação dos critérios de avanço

| Critério | Estado |
| --- | --- |
| Evidência externa consistente | **Não cumprido**: sinal positivo, mas IC inclui zero e sensibilidades não concordantes |
| Definição completamente congelada | Cumprido |
| Independência de um único doente | Cumprido na análise primária |
| Compatibilidade biológica do endpoint | Parcial, não plena |
| Evidência suficiente para dinâmica | Não cumprido |

Pelos critérios pré-especificados e pelos kill criteria do projeto, a Fase 2 não autoriza progressão para Level 2 external prediction, Digital Twin, simulação, classificador clínico ou virtual drug screening. O resultado negativo quanto ao avanço é retido como conclusão científica; P8 não será otimizado com ARTEMIS.

## Próximo teste cientificamente justificável

Uma futura tentativa exigiria uma coorte independente com endpoint mais homólogo e, para seleção versus plasticidade, amostras baseline/residual pareadas com DNA somático/CNV ortogonal. Essa aquisição constitui uma nova fase e deve ter protocolo próprio antes de abrir outcomes.
