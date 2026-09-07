# Definições congeladas dos programas

Versão: `v1`  
Rank: 9  
Seed: 11  
SHA-256 da definição: `c353b2b5edb47c18e960aeb10622dc673e6f6fdf24342deebd784a3982371aec`

Os scores são contínuos. Não existe threshold binário. Cada score é a média ponderada de `log1p(CPM)` pseudobulk, dividida pelo desvio-padrão congelado de cada gene. Os pesos somam um por programa.

| Programa | Rótulo estrutural | Genes principais | Coerência média | Correlação com total técnico | Redundância máxima |
|---|---|---|---:|---:|---:|
| P1 | housekeeping–citoesqueleto | SERF2, ACTG1, PPDPF, H3F3A, MYL6, TMSB4X, PFN1, NME2 | 0,177 | 0,74 | 0,68 |
| P2 | resposta imediata–stress | MALAT1, NEAT1, VMP1, SMG1, XIST, WSB1, RSRP1, SLC38A2 | 0,215 | 0,89 | 0,37 |
| P3 | proliferação–ciclo celular | H2AFZ, STMN1, TUBA1B, HMGN2, TUBB, UBE2C, TPI1 | 0,201 | 0,73 | 0,68 |
| P4 | immune-like–CD45 | PTPRC, IL7R, CD37, CD52, CD7, SRGN, LTB, RHOH | 0,196 | 0,53 | 0,26 |
| P5 | androgénico–secretor | HPX, NFATC4, TNNI3, LYPD3, MAGEC2, COL4A6, AR, KMO | 0,215 | 0,82 | 0,53 |
| P6 | apócrino–lipídico/metabólico | ALOX15B, ABCC11, SERHL2, ACSM1, DHRS2, MUCL1, ABCA12, KYNU | 0,261 | 0,57 | 0,42 |
| P7 | basal–mesenquimal/ECM | BGN, WIF1, SCRG1, PIK3R1, MEF2C, ECRG4, CRISPLD1, KRT14 | 0,180 | 0,82 | 0,46 |
| P8 | apresentação antigénica–interferão | B2M, CD74, HLA-B, HLA-C, HLA-A, IFITM3, HLA-DRA, UBD | 0,228 | 0,74 | 0,62 |
| P9 | luminal–secretor | TFF3, AGR3, CRISP3, CREB3L1, SCGB2A2, TSPAN1, AGR2, CLCA2 | 0,177 | 0,69 | 0,53 |

Coerência é a correlação de Pearson média entre os 30 genes com maior peso na amostra celular cega e balanceada. A correlação técnica é com a expressão total dos 2 000 genes selecionados por célula. A correlação elevada de P1, P2, P3, P5, P7 e P8 obriga a interpretar variação de score com cautela.

P2 e P3 são confundidores explícitos. P4 pode refletir contaminação, doublets residuais ou expressão immune-like genuína em células anotadas como malignas; não é tratado automaticamente como estado tumoral. Nenhum programa foi removido após inspeção.

Definições completas e pesos: `results/phase1/blind_discovery/program_definitions_v1.npz` e `program_definitions_v1.json`.
