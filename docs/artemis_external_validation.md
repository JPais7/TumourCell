# Validação externa de P8 em ARTEMIS

Data: 2026-09-06  
Estado: análise executada depois da congelação do protocolo  
Conclusão: **INCONCLUSIVO**

## Pergunta e interpretação do endpoint

Foi aplicada sem retuning a definição P8 congelada em GSE246613. O contraste ARTEMIS é P8 no baseline em doentes que atingiram pCR menos P8 no baseline em doentes com doença residual (RD). A direção prevista foi positiva.

Este é um teste externo parcial de associação com sensibilidade terapêutica. Não é uma replicação literal do contraste longitudinal R2–NR após pembrolizumab da discovery: ARTEMIS contém tecido molecular pré-tratamento e outro regime terapêutico.

## Fronteira anti-leakage

Antes de abrir o outcome foram congelados:

- protocolo: SHA-256 `583e704151fc3f17536e9b9b9912f3922391aab9ef4d6a3afbbded3087350897`;
- definição P8 exportada: SHA-256 `185cde59ee15f1fac893b8367cfd1320dd6e1559d8e4894ab060cebfc4ed744e`;
- definição original da Fase 1: SHA-256 `c353b2b5edb47c18e960aeb10622dc673e6f6fdf24342deebd784a3982371aec`.

Não foi reexecutado NMF, não foram escolhidos genes, pesos, rank, seed, normalização ou thresholds com base em ARTEMIS.

## Auditoria do dataset

| Campo | Observado |
| --- | --- |
| Objeto | CELLxGENE `e94bd3cc-6271-424a-baac-12f8eb320a0e` |
| Ficheiro | H5AD, 4 522 982 499 bytes; SHA-256 `d40a2b210ddf8fe0bd46118900412493c6d11d3429ffae4c311f7469f7d10da0` |
| Células | 427 823 |
| Doentes | 101 |
| Células tumorais | 49 275, anotação final `author_cell_type=Tumor` / `cell_type=abnormal cell` |
| Doentes com células tumorais | 97 |
| Outcome entre estes doentes | 45 pCR, 34 RD, 18 `Excluded` |
| Plataforma | 10x 3′ v3: 405 118 células; 10x 3′ v2: 22 705 |
| Momento molecular | Pré-tratamento |
| Matriz usada | `raw/X`, counts não negativos com valores inteiros representados em `float32` |
| Genes | 32 354 símbolos únicos |
| Amostras por doente | Uma identidade de dador no objeto; counts tumorais agregados por dador |

A anotação tumoral é derivada da estratégia publicada de aneuploidia baseada em RNA. É utilizada para inclusão celular, não como prova de clone.

## Compatibilidade de P8

- 1 927/2 000 genes presentes (requisito ≥1 600);
- 98,46% do peso total presente (requisito ≥80%);
- 49/50 genes principais presentes (requisito ≥40); apenas `PLAAT4` ausente;
- os pesos presentes foram renormalizados, conforme pré-especificado;
- cobertura: **APROVADA**.

## Análise primária

Foram somados counts tumorais por doente, calculado `log1p(CPM)`, aplicada a escala génica e os pesos P8 congelados, e padronizado o score nos doentes tecnicamente elegíveis antes da estratificação por outcome.

| Estimativa | Resultado |
| --- | ---: |
| Doentes analisados | 64 (36 pCR; 28 RD) |
| Células tumorais mínimas | 50 por doente |
| Média pCR | 0,259 SD ARTEMIS |
| Média RD | 0,072 SD ARTEMIS |
| Diferença pCR−RD | **+0,187 SD** |
| IC95% bootstrap | **−0,242 a +0,607** |
| Diferença de medianas | +0,541 SD |
| Hedges g | +0,208 |

**Observado:** o estimador pontual tem a direção prevista, mas o intervalo inclui efeitos negativos, zero e efeitos positivos moderados. Logo, não satisfaz o critério congelado de replicação parcial.

## Influência por doente

- A direção foi positiva em 100% das 64 análises leave-one-patient-out.
- Intervalo dos efeitos leave-one-out: +0,143 a +0,269 SD.
- Mais influentes: P81 (pCR; efeito sem o doente +0,269), P55 (RD; +0,242) e P99 (RD; +0,241).

**Inferido:** a análise primária não depende de um único doente. Isto não compensa a imprecisão nem a instabilidade entre especificações.

## Sensibilidades pré-especificadas

| Especificação | n pCR/RD | Efeito pCR−RD (SD) | IC95% | LOO positivo |
| --- | ---: | ---: | ---: | ---: |
| ≥25 células, log-CPM | 38/31 | +0,282 | −0,125 a +0,673 | 100% |
| **≥50 células, log-CPM (primária)** | **36/28** | **+0,187** | **−0,242 a +0,607** | **100%** |
| ≥100 células, log-CPM | 30/25 | −0,010 | −0,478 a +0,436 | 38,2% |
| ≥50 células, log-TP10k | 36/28 | +0,305 | −0,141 a +0,742 | 100% |
| ≥50 células, sqrt-CPM | 36/28 | +0,268 | −0,178 a +0,705 | 100% |
| ≥50, ajuste conjunto por log-library, nº células e versão 10x | 36/28 | −0,002 | não usado para decisão primária | — |

No conjunto primário, seis doentes foram 10x v2 (2 pCR, 4 RD) e 58 foram v3 (34 pCR, 24 RD). O score correlacionou-se com log-library (r=−0,43) e número de células tumorais (r=−0,40); log-library e número de células correlacionaram-se entre si (r=0,75). O ajuste conjunto deve ser lido como sensibilidade conservadora, com colinearidade relevante, não como nova análise primária.

**Observado:** a direção é preservada nas normalizações e no threshold de 25 células, mas desaparece no threshold de 100 e no ajuste técnico. Todos os IC95% não ajustados incluem zero.

## Comparação com a discovery

Na discovery, o efeito P8 após pembrolizumab foi +0,343 SD (IC95% −0,078 a +0,781) e era vulnerável à remoção de Patient45. Em ARTEMIS, o estimador primário é menor (+0,187 SD), na mesma direção, e não inverte no leave-one-out. As magnitudes e intervalos são compatíveis em sentido amplo, mas os endpoints não são homólogos e ambas as estimativas são imprecisas.

## Veredicto congelado

**INCONCLUSIVO.**

- Não é **REPLICAÇÃO PARCIAL**, porque o IC95% primário inclui zero e a direção não se conserva em todas as sensibilidades principais.
- Não é **NÃO REPLICAÇÃO**, porque a estimativa primária tem a direção prevista e o IC95% não está inteiramente abaixo de zero.
- O resultado é compatível com H1 (associação biológica real), H5 (efeito técnico) e H6 (efeito fraco/específico da coorte); não as distingue.
- Não há suporte para transformar P8 num biomarcador clínico nem para inferir mecanismo, seleção ou plasticidade.

## Artefactos reproduzíveis

- Script: `src/validate_p8_artemis.py`
- Resultado completo: `results/phase2/artemis_p8_validation.json`
- Scores paciente-nível: `results/phase2/artemis_p8_patient_scores.npz`
- Figura: `results/figures/phase2_p8_artemis.png`
- Seed bootstrap/jitter: 20260906; 10 000 reamostragens estratificadas por outcome.

## Proveniência pública

- Artigo: https://www.nature.com/articles/s41586-026-10469-9
- Repositório e código dos autores: https://github.com/navinlabcode/tnbc-chemo
- CELLxGENE: https://cellxgene.cziscience.com/e/6f9de485-58cd-4342-bfc4-b3d3dd223aa8.cxg/

