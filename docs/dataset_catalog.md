# Catálogo sistemático de datasets — Atlas de Estados Malignos TNBC

Data da pesquisa: 2026-09-06  
Âmbito: tecido tumoral humano, sc/snRNA-seq transcriptoma-wide, TNBC, com counts e identidade de doente  
Fontes pesquisadas: GEO/SRA, CELLxGENE Discover, artigos originais e repositórios dos autores

## Método de levantamento

Foram combinadas pesquisas por `TNBC`, `triple-negative breast cancer`, `single-cell`, `single-nucleus`, `scRNA-seq`, `malignant`, `atlas`, `neoadjuvant` e `residual`, seguidas de reconciliação por accession, publicação e doente. Reanálises dos mesmos accessions são registadas, mas não contam como coortes independentes. Os critérios de elegibilidade foram congelados em `docs/state_definition.md` antes de qualquer novo download da Fase 3.

O catálogo é sistemático dentro das fontes públicas indexadas e da data de corte, mas não é uma revisão PRISMA exaustiva de dados controlados ou ainda não publicados.

## Coortes diretamente avaliadas para o atlas primário

| Estudo / accession | Doentes e células | Malignas | Plataforma | Tratamento/momento | Counts e metadados | CNV / spatial / follow-up | Papel e decisão pré-análise |
| --- | --- | ---: | --- | --- | --- | --- | --- |
| Shiao, GSE246613 | 50 inscritos; 33 com série scRNA completa; 171 746 células não imunes no objeto usado | 50 185 cancer cells | 10x scRNA/snRNA | baseline; pembrolizumab; pembrolizumab+RT | H5AD raw e mapas doente/tempo locais | CNV não fornecido; spatial protein; resposta longitudinal | **Construção**; programas Fase 1 congelados, sem repetir NMF |
| Yan/ARTEMIS, PRJNA1041570; CELLxGENE `e94bd3cc…` | 101 doentes; 427 823 células | 49 275 | 10x 3′ v2/v3 | pré-tratamento, NAC AC-T | H5AD raw validado; doente e anotação tumoral | RNA-CNV/CopyKAT; Xenium/Visium; pCR/RD | **Construção** outcome-blind; passa dimensão/acesso |
| Wu, GSE176078 | 26 tumores: 11 ER+, 5 HER2+, 10 TNBC; ~100 mil células totais | 24 489 epiteliais/cancer no depósito CZI, subconjunto TNBC a auditar | 10x 3′ v2 e 5′ v1 | primário não tratado | GEO e H5AD CZI; IDs de dador/subtipo | inferred CNV; Visium; sem follow-up terapêutico | **Construção**, apenas TNBC; elegibilidade final após auditoria ≥3×50 malignas |
| Karaayvaz, GSE118389 | 6 TNBC; 1 534 células | maioria tumoral; número final após QC | Fluidigm C1 / full-length RNA-seq | cirurgia primária, sem tratamento sistémico reportado antes da colheita | GEO: counts RSEM 9 MB, TPM 11,5 MB; IDs por célula/doente | RNA-CNV e WES na SuperSeries; sem spatial/follow-up | **Construção**; contraste de plataforma, condicionado a ≥500 malignas e QC |
| Pal, GSE161529 | 52 dadores totais, 32 tumores; 421 761 células reportadas; subconjunto TNBC a auditar | epiteliais malignas no atlas | 10x | maioritariamente treatment-naive; primário/nodo em subset | GEO/CELLxGENE; doente, subtipo e tecido | CNV não primário; sem resposta; nós pareados em subset | **Validação reservada**; nunca usada para limiares/assinaturas |

## Coortes candidatas ou contextuais

| Estudo | Informação útil | Avaliação pelos critérios congelados |
| --- | --- | --- |
| Bassez 2021, EGAS00001004809 / EGAD00001006608 | 40 doentes de vários subtipos, pre/on anti-PD1, scRNA/CITE/TCR/exome | **NÃO DETERMINÁVEL para inclusão pública imediata**: acesso controlado e subconjunto TNBC/maligno necessita autorização e auditoria |
| GSE266919 | TNBC avançado sob taxanos ± atezolizumab; integração de dados novos e publicados | Candidato de sensibilidade futura; risco de duplicação e mistura de fontes. Não substitui coorte pré-atribuída nesta versão |
| GSE299631 / BRCA1 residual 2026 | modelos BRCA1 e componente humano espacial/single-cell | **NÃO DETERMINÁVEL** até separar espécie, doente e counts humanos; não usado nesta versão |
| GSE276609 / PDX platinum | TNBC PDX com scRNA e single-cell CN longitudinal | **EXCLUÍDO do atlas humano primário**; contexto mecanístico futuro |
| GSE169246 | 489 490 células CD45+ em TNBC tratada | **EXCLUÍDO**: não contém compartimento maligno |
| GSE205472 | HR+ NAC, 8 doentes scRNA | **EXCLUÍDO**: subtipo não-TNBC; controlo futuro de especificidade |
| GSE114727 / FELINE | ER+/HER2− serial sob letrozole±ribociclib | **EXCLUÍDO**: subtipo não-TNBC |
| Wang 2023, Hammerl 2021, Bassiouni 2023 | spatial/imaging ou bulk com contexto TNBC | **EXCLUÍDOS da definição de estados**: sem scRNA maligno transcriptoma-wide independente adequado; podem localizar estados posteriormente |
| Kim 2018, PRJNA396019 | 6 862 sc/snRNA, 8 casos detalhados; 900 scDNA; serial | Cientificamente valioso, mas raw público ~1,4 TiB e matriz/processamento ainda não resolvidos; **NÃO DETERMINÁVEL nesta versão** em vez de baixar limiares de auditabilidade |

## Diferenças técnicas antecipadas

- ARTEMIS e Wu/Pal usam tecnologias 10x mas versões distintas; Karaayvaz usa full-length/Fluidigm C1.
- GSE246613 inclui exposição terapêutica e momentos repetidos; ARTEMIS, Wu, Karaayvaz e o núcleo da validação Pal são predominantemente baseline.
- Identificação maligna varia entre anotação publicada, inferência RNA-CNV e compartimento epitelial; será conservada por coorte e auditada, não harmonizada por outcome.
- Profundidade por célula e número de células por doente variam por ordens de grandeza; o subsampling por doente é obrigatório.
- A mesma matriz integrada que contenha células Wu e Pal pode ser usada apenas como veículo de acesso aos counts; as fontes são separadas por prefixo de dador e Pal permanece oculto até congelação do atlas.

## Estado do acesso na data de corte

| Dataset | Estado |
| --- | --- |
| GSE246613 | local, checksum verificado |
| ARTEMIS | local, checksum `d40a2b210ddf8fe0bd46118900412493c6d11d3429ffae4c311f7469f7d10da0` |
| GSE118389 | ficheiros públicos pequenos identificados; ainda não descarregados na fronteira de congelação |
| GSE176078 | URL H5AD público identificado; ainda não descarregado na fronteira de congelação |
| GSE161529 | holdout identificado; não abrir/analisar antes do freeze das assinaturas |

## Fontes primárias

- ARTEMIS: https://www.nature.com/articles/s41586-026-10469-9 e https://cellxgene.cziscience.com/collections/0a117356-ecaa-4b82-a454-c15a4c9ec507
- GSE246613: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246613
- Wu/GSE176078: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078 e https://cellxgene.cziscience.com/collections/65db5560-7aeb-4c66-b150-5bd914480eb8
- Karaayvaz/GSE118389: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118389
- Pal/GSE161529: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE161529
- Atlas epitelial integrado usado apenas para localizar ativos: https://cellxgene.cziscience.com/collections/9432ae97-4803-4b9f-8f64-2b41e42ad3cb

