# Protocolo congelado de validação externa — Fase 2

Data de congelação: 2026-09-06  
Estado: escrito antes do download do objeto ARTEMIS e antes de abrir o outcome  
Coorte de descoberta: GSE246613  
Coorte externa: Yan/ARTEMIS 2026

## Definição imutável de P8

- Programa: P8, apresentação antigénica/interferão.
- Definição de origem: `program_definitions_v1.npz`.
- SHA-256: `c353b2b5edb47c18e960aeb10622dc673e6f6fdf24342deebd784a3982371aec`.
- Rank/seed de origem: 9/11.
- Genes, pesos e escalas completos: `results/phase2/p8_definition_frozen.json`.
- Genes de maior peso: B2M, CD74, HLA-B, HLA-C, HLA-A, IFITM3, HLA-DRA, ACTG1, UBD, TAGLN2, S100A6, BST2, TPT1, IFI27, PPIA, MMP7, CFB, HLA-DRB1, ANXA2, RARRES1, SAA1, LY6E, IFITM1, KRT19, UBC, PSMB9, PLAAT4, WARS, S100A9, PRDX1, MGST1, CSTB, ALDOA, PFN1, PPIB, CFL1, UBE2L6, TPI1, CALR, IFI6, S100A14, TYMP, MYL6, SELENOM, PDZK1IP1, HLA-DRB5, ASS1, CXCL10, S100A16 e HSP90AB1.

Não haverá NMF, seleção de genes, escolha de outro programa ou alteração de pesos em ARTEMIS.

## Compatibilidade e endpoint

ARTEMIS é pré-tratamento e avalia resposta à quimioterapia neoadjuvante. O endpoint primário é pCR versus doença residual (RD). Isto é uma validação **parcial** da associação baseline de P8 com sensibilidade terapêutica, não uma replicação literal da alteração R2–NR após radioimunoterapia.

Contraste primário congelado: média de P8 em pCR menos média em RD. Direção prevista: positiva.

## Unidade e população

- Unidade inferencial: doente.
- Apenas células malignas/aneuploides segundo a anotação final publicada pelos autores.
- Pseudobulk de counts por doente antes de calcular o score.
- Uma observação por doente; amostras técnicas do mesmo doente são agregadas.
- Inclusão primária: pelo menos 50 células malignas e endpoint pCR/RD conhecido.
- Sensibilidades: pelo menos 25 e pelo menos 100 células malignas.
- Não excluir doentes com base no valor de P8 ou no outcome.

## Scoring e genes ausentes

1. Somar counts brutos por doente.
2. Calcular `log1p(counts/library_size × 1 000 000)`.
3. Dividir cada gene pela escala congelada da descoberta.
4. Aplicar os pesos P8 congelados e somar.
5. Padronizar o score entre todos os doentes tecnicamente elegíveis, antes de estratificar por outcome.

Genes ausentes recebem contribuição zero. Os pesos só são renormalizados sobre genes presentes se estiverem disponíveis simultaneamente:

- pelo menos 80% do peso total P8;
- pelo menos 40 dos 50 genes de maior peso;
- pelo menos 1 600 dos 2 000 genes da definição.

Se algum requisito falhar: **NÃO DETERMINÁVEL**.

## Modelo estatístico

- Estimador primário: diferença de médias pCR−RD em unidades de desvio-padrão ARTEMIS.
- Incerteza: bootstrap não paramétrico de 10 000 reamostragens de doentes dentro de cada grupo, seed 20260906.
- Reportar também diferença de medianas e Hedges g.
- Leave-one-patient-out: direção, variação do efeito e três doentes mais influentes.
- Análise primária não ajustada.
- Sensibilidades: ajuste linear por `log(library_size)`, número de células malignas e batch/coorte quando identificável; log-TP10k; sqrt-CPM; thresholds 25/100.
- Missing endpoint: exclusão explícita com contagem; sem imputação.
- Batch completamente confundido com outcome implica **NÃO DETERMINÁVEL** para inferência ajustada.

## Critérios de interpretação

### REPLICAÇÃO PARCIAL

Efeito pCR−RD positivo, IC95% totalmente acima de zero, cobertura genética/QC aprovados e direção preservada em pelo menos 90% das análises leave-one-out e nas sensibilidades principais.

### NÃO REPLICAÇÃO

Efeito negativo com IC95% totalmente abaixo de zero e direção negativa estável nas sensibilidades; ou incompatibilidade consistente numa amostra com precisão suficiente.

### NÃO DETERMINÁVEL / INCONCLUSIVO

IC95% inclui zero, endpoint apenas parcialmente comparável, cobertura insuficiente, confounding técnico não resolúvel, amostra efetiva pequena ou dependência material de poucos doentes.

Não serão criados thresholds pós-hoc. Uma associação do metaprograma interferão definido pelos autores não substitui o teste dos pesos P8 congelados.

## Proveniência técnica

- Artigo: https://www.nature.com/articles/s41586-026-10469-9
- Repositório: https://github.com/navinlabcode/tnbc-chemo
- Commit auditado: `fbe1dd3b1db05dd5e9f02485a45fdd438d6af9fb`
- CELLxGENE: https://cellxgene.cziscience.com/e/6f9de485-58cd-4342-bfc4-b3d3dd223aa8.cxg/

Este documento e `p8_definition_frozen.json` constituem a fronteira anti-leakage da Fase 2.
