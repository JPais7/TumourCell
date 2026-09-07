# Seleção versus plasticidade

## Resultado

**INDETERMINADO — os dados atuais não permitem distinguir seleção de reprogramação transcricional.**

## O que é observado

Todos os nove programas têm expressão mensurável nos pseudobulks baseline. Após a adição de radioterapia, P3 diminui em R2 e NR; P2 e P4 diminuem em NR. P8 apresenta um contraste R2–NR positivo, mas com incerteza ampla e intervalo compatível com ausência de diferença.

## Porque isto não identifica seleção

Um score pseudobulk maior pode resultar de alteração da abundância relativa de estados malignos preexistentes, alteração de expressão dentro das células, qualidade/profundidade da biopsia ou composição residual. Não foram estimadas proporções de estados com uma regra de classificação congelada e não existe lineage tracing entre biopsias.

## Porque isto não identifica plasticidade

O H5AD não contém CNV, clones ou genómica longitudinal. `obsm` contém apenas `X_scVI` e `X_umap`; não existe uma variável que permita ligar células baseline e pós-tratamento à mesma população clonal. Um aumento transcricional após tratamento não demonstra reprogramação dentro de clones persistentes.

## Classificação da evidência

- **OBSERVADO:** alterações de scores pseudobulk entre biopsias do mesmo doente.
- **INFERIDO:** os programas resumem eixos de covariação da expressão maligna.
- **INCERTO:** se P8 acompanha especificamente resposta tardia à radioterapia.
- **NÃO IDENTIFICÁVEL:** seleção versus plasticidade; expansão versus indução; persistência clonal.

O próximo teste mínimo para esta pergunta requer CNV inferido validado por doente/tempo ou dados genómicos emparelhados. Mesmo CNV inferido permitiria apenas “compatível com alteração dentro de estrutura copy-number persistente”, não lineage tracing direto.
