# Fase 4 — Gate para digital twin molecular

Data: 8 de setembro de 2026.

## Veredicto

**Não está autorizada uma simulação individual de tratamentos ou perturbação virtual.** Está autorizada uma representação molecular observacional do doente, com incerteza explícita.

Este veredicto aplica os critérios de paragem já definidos no projeto: não avançar para perturbação virtual enquanto a previsão externa de nível 2 não generalizar.

## Teste externo

Foram avaliadas as duas assinaturas prioritárias nas 42 amostras basais independentes do GSE260693, sem refitting. O resultado favorável foi ausência de tumor residual.

| Assinatura | AUROC | IC bootstrap 95% | Gate: limite inferior >0,50 |
|---|---:|---:|---|
| State_03 | 0,613 | 0,430–0,786 | Falha |
| P8 | 0,647 | 0,471–0,817 | Falha |

Ambas têm direção promissora, mas a incerteza permite desempenho equivalente ao acaso. A associação combinada entre coortes não substitui validação preditiva externa.

## Estado dos gates

1. Transporte da representação congelada: **parcial** — cobertura génica elevada, mas bulk não é específico das células malignas.
2. Previsão externa ao nível do doente: **falha** — intervalos AUROC atravessam 0,50.
3. Transições condicionadas pelo tratamento: **bloqueado** — braço e pCR do NeoTRIP indisponíveis; regimes heterogéneos no GSE260693.
4. Dinâmica clonal: **bloqueado** — não existe validação longitudinal clone-resolved adequada.
5. Validação espacial: **bloqueado** — falta coorte TNBC whole-transcriptome com outcome e dimensão suficiente.

## Artefacto permitido

A ficha molecular v0.1 contém:

- valores de State_01–State_04 e P8;
- percentis dentro da mesma coorte e momento;
- contexto terapêutico, quando disponível;
- incerteza de medição, domain shift e ambiguidade biológica;
- modalidades ausentes;
- restrição explícita a descrição observacional e investigação de coortes.

Não contém probabilidade clínica, tratamento recomendado, transição causal ou sobrevivência simulada.

## Critério para desbloquear o simulador

É necessário obter uma coorte externa adicional ou os metadados clínicos do NeoTRIP e demonstrar, com análise congelada ao nível do doente, discriminação cujo IC95% fique acima de 0,50, acompanhada de calibração adequada. A seguir será necessário estimar transições separadamente por tratamento e validar essas transições numa coorte independente.
