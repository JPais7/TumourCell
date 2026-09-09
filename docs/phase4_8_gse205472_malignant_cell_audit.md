# Phase 4.8 — Observabilidade de células malignas no GSE205472

Foi inspeccionado o registo GEO/SOFT e os links públicos para matriz, features e
barcodes. O dataset é scRNA-seq de biopsias tumorais, mas não foi encontrada uma
annotation cell-level `malignant`, nem CNV inferida, clone, mutação ou ground
truth patológico celular nos artefactos inspeccionados.

`EPCAM`/queratinas, exclusão imune/estromal e o rótulo “tumor sample” não são
aceites como malignidade. Outcome e tratamento não entram na regra. Assim,
`epithelial_cell`, `tumour_like_cell` e `malignant_cell` permanecem conceitos
distintos; a classificação actual é `UNRESOLVED` para as células e
`MALIGNANT_CELL_OBSERVABLE_NOT_VERIFIED` para o dataset.

Há cinco pares pré/pós com patient IDs, mas a estabilidade de malignidade
pré/pós não é avaliável sem uma regra cell-level congelada. A Phase 5 continua
bloqueada; não foram treinadas dinâmicas nem feitas recomendações clínicas.
