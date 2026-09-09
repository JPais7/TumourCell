# Phase 4.9 — Observabilidade global de células malignas

Foram auditados GSE205472, GSE176078, GSE246613, SRP114962, NeoLetExe e
ARTEMIS. A matriz distingue ground truth cell-level, inferência reprodutível,
evidência parcial e ausência de observabilidade. EPCAM/KRT, amostra tumoral,
exclusão imune/estromal e CNV inferida não são automaticamente malignidade.

Não foi encontrado dataset com ground truth E4/E5 suficiente para desbloquear a
via de dinâmica. Alguns recursos podem suportar inferência futura (GSE176078,
SRP114962, NeoLetExe e ARTEMIS), mas dependem de annotations/CNV e validação
independente. GSE205472 permanece sem observabilidade maligna cell-level; a sua
longitudinalidade discreta não resolve essa lacuna.

Classificação global: **PUBLIC_MALIGNANT_OBSERVABILITY_PARTIAL**. Outcome e
tratamento não entram na definição. A Phase 5 permanece bloqueada; não foram
implementados modelos temporais, simulação ou recomendações clínicas.

## Annotation semantics

Uma annotation de autor só é evidência de malignidade quando o tipo é
explicitamente `malignant`, a proveniência é E4/E5 e existe validação adequada.
`epithelial`, `tumour_like`, `tumour_sample`, EPCAM/KRT e exclusão imune são
apenas evidência parcial. CNV inferida só pode atingir inferência reprodutível
quando o pipeline está explicitamente suportado; nunca é ground truth automática.
Outcome e tratamento não podem promover uma classificação.
