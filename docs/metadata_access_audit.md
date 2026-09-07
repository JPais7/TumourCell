# Auditoria de metadados e acesso — Gate 1

Data: 2026-09-06

## Resultado

O Gate 1 fica **parcialmente aprovado**. ARTEMIS/Yan e GSE246613 têm dados processados públicos e estrutura suficiente para funções específicas. O estudo de Kim tem dados brutos públicos, mas a via de reprocessamento integral é impraticável para a fase inicial sem uma matriz processada e um mapa célula–doente–tempo. O GSE205472 é utilizável como controlo exploratório, mas não constitui oito pares completos.

## Kim 2018 — dinâmica primária de quimioterapia

- SRA confirmado: `SRP114962`, resolvido pelo ENA para `PRJNA396019`.
- O ENA devolve 6 633 runs/títulos de amostra e aproximadamente **1 386 GiB** de FASTQ comprimido.
- O artigo descreve 6 862 núcleos RNA de oito doentes, com RNA longitudinal detalhado em seis doentes de persistência e extinção, e scDNA em oito doentes.
- Os identificadores codificam doente e tempo (`KTN…0`, `…2`, `…OP`), mas o SRA mistura snRNA, scDNA, exoma, sangue e amplicões.
- Não foi confirmada nesta auditoria uma matriz de expressão processada oficial, diretamente descarregável, acompanhada pelo dicionário completo de células.

**Decisão:** não descarregar 1,4 TiB. Procurar uma matriz processada nos suplementos/reanálises reproduzíveis ou contactar os autores. Se não existir, gerar apenas uma lista filtrada de runs de snRNA e estimar custo antes de qualquer download. Isto é um bloqueio operacional, não uma falsificação científica.

Fontes: [artigo](https://pubmed.ncbi.nlm.nih.gov/29681456/), [SRA](https://www.ncbi.nlm.nih.gov/sra/?term=SRP114962), [ENA](https://www.ebi.ac.uk/ena/browser/view/PRJNA396019).

## GSE246613 — validação longitudinal por mudança de regime

- 50 doentes foram incluídos no ensaio, mas o depósito contém dados moleculares humanos de **34 doentes**.
- O SOFT contém 252 registos humanos: 104 scRNA-seq, 46 snRNA-seq e 102 scTCR-seq; existem ainda 14 registos murinos, que ficam excluídos da análise clínica.
- 33 dos 34 doentes têm scRNA-seq nos três momentos A/B/C: baseline, após pembrolizumab e após pembrolizumab+radioterapia. Um doente tem dois momentos.
- 15 doentes têm snRNA nos três momentos; 31 têm scTCR nos três momentos.
- Os ficheiros processados principais estão públicos: 1,0 GB (imunes), 1,4 GB (não imunes), 376,6 MB (objeto combinado) e 14,9 MB (TCR).
- O SOFT não inclui a classe de resposta. Essa variável terá de ser recuperada do objeto H5AD, suplemento do artigo ou fonte clínica antes de abrir os dados de expressão para o teste bloqueado.

**Decisão:** conjunto forte para falsificação longitudinal e especificidade de regime, mas não substitui a descoberta de quimioterapia. Autorizar apenas a transferência do H5AD mínimo necessário depois de confirmar onde reside a resposta.

Fontes: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246613), [artigo](https://doi.org/10.1016/j.ccell.2023.12.012).

## GSE205472 — controlo HR+ de quimioterapia

- 13 amostras scRNA: 8 pré-NAC e 5 pós-NAC.
- Oito doentes no baseline; apenas P01, P02, P04, P05 e P07 têm amostra pós-tratamento depositada.
- Resposta no baseline: 3 sensíveis e 5 resistentes. Entre os cinco pares: 2 sensíveis e 3 resistentes.
- A matriz processada contém **91 235 barcodes/células**.
- O depósito fornece uma matriz 10x combinada, mas o mapeamento dos prefixos de barcode para as 13 amostras e as anotações celulares não está documentado de forma suficiente no SOFT.
- O registo GEO continua sem citação de publicação associada e indica contacto com o investigador principal para condições de utilização.

**Decisão:** manter como controlo negativo/exploratório, não como validação central. Exigir mapa de barcode e confirmar os pares ausentes antes de análise.

Fonte: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE205472).

## GSE169246 — dinâmica imunitária sob paclitaxel ± atezolizumab

- O depósito contém 489 490 barcodes de scRNA-seq em 78 amostras: 48 de sangue e 30 de tumor, provenientes de 22 doentes.
- Existem 16 doentes com pelo menos uma amostra tumoral e 12 com tumor emparelhado pré/pós: seis tratados com paclitaxel e seis com paclitaxel+atezolizumab. Apenas um tem tumor nos três momentos pré/pós/progressão.
- O desenho e o depósito especificam enriquecimento CD45+. Não existe um compartimento maligno adequado para testar seleção, plasticidade ou programas residuais de células tumorais.
- O GEO fornece matriz processada completa: barcodes (2,0 MB), features (98,3 KB) e counts (1,9 GB). O manifesto foi construído apenas a partir do SOFT e dos barcodes; os counts não foram descarregados nem lidos.
- A resposta RECIST foi localizada num recurso secundário derivado do artigo para 21 dos 22 doentes: 9 R e 12 NR; P028 não tem resposta nesse recurso. Este mapa não será tratado como metadado primário até ser reconciliado com o suplemento original.

**Decisão:** rejeitar como substituto de Kim para dinâmica maligna. Manter como validação mecanística imunitária alinhada com paclitaxel, com análises estratificadas por braço e compartimento.

Fontes: [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169246), [artigo](https://doi.org/10.1016/j.ccell.2021.09.010), [CellResDB](https://cellknowledge.com.cn/cellresponse/php_mysql/datasetDetail.php?id=GSE169246_TNBC_1).

## Yan/ARTEMIS 2026 — validação de resposta e espaço

- 101 doentes com scRNA-seq; 427 857 células totais; 97 doentes com células malignas e 49 275 células malignas.
- 44 doentes com Xenium; 39 com scRNA e espaço emparelhados.
- Resposta disponível para 89 doentes: 50 pCR e 39 doença residual; 19 sem endpoint.
- Dados públicos em `PRJNA1041570`, CELLxGENE e código/suplementos no repositório dos autores.
- É um coorte exclusivamente pré-tratamento para o ensaio molecular, pelo que valida associação/predição e contexto espacial, não emergência temporal.

**Decisão:** manter completamente bloqueado para validação. Não usar para definir programas, genes, sinais ou thresholds.

Fontes: [Nature](https://www.nature.com/articles/s41586-026-10469-9), [código](https://github.com/navinlabcode/tnbc-chemo).

## Decisão operacional seguinte

1. Recuperar ou reconstruir o dicionário de células do estudo Kim sem descarregar todos os FASTQ. O TISCH2 confirma uma versão processada de 2 472 células/8 doentes, mas o link de download exposto pelo site devolvia 404 em 2026-09-06.
2. A variável de resposta do GSE246613 foi localizada no H5AD não imune e foi produzido um manifesto clínico separado da expressão.
3. Obter o mapa barcode–amostra do GSE205472; se não for recuperável, substituir este controlo por um coorte HR+ mais reprodutível.
4. Só depois congelar o protocolo estatístico e iniciar qualquer análise de expressão.

## Critério de pivot imediato

Se não for possível obter uma matriz Kim processada com metadados suficientes, o projeto não deve fingir uma análise clone-aware. O pivot correto é uma análise de redistribuição longitudinal de programas em GSE246613, com a questão explicitamente reformulada para resposta a pembrolizumab+radioterapia, ou a procura de outro coorte de quimioterapia TNBC com pares públicos.

## Atualização GSE246613 após inspeção do H5AD

- O objeto combinado de 376,6 MB contém apenas dados murinos e não deve ser usado para inferência clínica.
- O objeto humano não imune tem 171 746 células, incluindo 50 185 anotadas como malignas.
- Resposta por doente após colapsar os dois tumores do doente 03: 9 R1, 14 R2 e 11 NR.
- Doentes com células malignas nos três momentos: 31 com pelo menos 1; 21 com pelo menos 10; 14 com pelo menos 25; 12 com pelo menos 50; 8 com pelo menos 100.
- A elegibilidade é, portanto, sensível ao limiar celular e tem de ser definida antes de calcular qualquer score biológico.
