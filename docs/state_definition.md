# Definição e protocolo congelado de estados malignos — Fase 3

Data de congelação inicial: 2026-09-06  
Estado: definido antes de descarregar novas coortes para a Fase 3 e sem consultar outcomes para definir estados

## Objetivo

Testar se existe um vocabulário mínimo de configurações transcricionais malignas recorrentes em coortes humanas independentes de TNBC. A hipótese nula é que as configurações são específicas de doente, coorte, tratamento ou plataforma.

Esta fase não procura resposta clínica, não reabre P8 e não altera qualquer artefacto das Fases 1–2.

## Critérios de elegibilidade de datasets

Uma coorte é elegível para descoberta/replicação de estados se cumprir todos os critérios:

1. tecido tumoral humano invasivo com pelo menos três doentes TNBC identificáveis;
2. scRNA-seq ou snRNA-seq transcriptoma-wide com matriz de counts pública;
3. identificador de doente e de amostra recuperável;
4. células malignas anotadas pelos autores ou identificáveis por uma regra publicada independente de outcomes;
5. pelo menos 500 células malignas no total e pelo menos três doentes com ≥50 células malignas;
6. pelo menos 5 000 genes medidos e símbolos/identificadores convertíveis de forma não ambígua;
7. metadados suficientes para plataforma, tecido e momento de colheita;
8. licença/acesso que permita reanálise auditável.

Critérios de exclusão:

- linhas celulares, organoides, PDX ou modelos murinos sem tecido humano primário analisável separadamente;
- coortes exclusivamente imunes, bulk ou spatial sem resolução celular transcriptoma-wide;
- estudos sem counts, sem doente ou sem identificação maligna defensável;
- reanálises/integrações de dados já incluídos, que não contam como replicação independente;
- amostras metastáticas misturadas com primárias sem metadados que permitam estratificação;
- coortes com menos de três doentes elegíveis após QC.

Falhar um critério gera `EXCLUÍDO` ou `NÃO DETERMINÁVEL`, nunca flexibilização pós-hoc.

## Papéis pré-atribuídos

- **Construção:** GSE246613, ARTEMIS, GSE176078/Wu e GSE118389/Karaayvaz, se passarem a auditoria técnica.
- **Validação completamente reservada:** GSE161529/Pal. Não participa na escolha de rank, matching, limiares ou assinaturas.
- **Não elegíveis para o atlas maligno primário:** GSE169246 (CD45+), datasets exclusivamente spatial/imaging, PDX e coortes não-TNBC.
- Outcomes, sobrevivência e resposta são ocultados durante definição e matching. Tratamento/plataforma são usados apenas depois para caracterizar dependência.

Se Wu ou Karaayvaz falharem acesso/QC, não serão substituídos com base no efeito observado; qualquer substituição exige uma versão datada deste protocolo antes da análise da nova coorte.

## Definição operacional de estado

Um candidato intra-coorte é um programa não negativo de coexpressão maligna, estimado separadamente dentro de cada coorte. Um **estado recorrente** é uma família de candidatos correspondentes em pelo menos duas coortes de construção e que passa estabilidade intra e inter-coorte. `State_01`, `State_02`, … são identificadores neutros atribuídos só após o matching.

Estados descrevem configurações contínuas de expressão. Hard clustering de células não é a definição primária; frequência celular é uma caracterização secundária baseada em uso relativo.

## Pré-processamento por coorte

1. usar counts brutos e apenas células malignas segundo a anotação publicada;
2. excluir genes mitocondriais e ribossomais da seleção de features, mantendo-os para QC;
3. exigir ≥500 UMIs e ≥300 genes por célula quando as métricas estão disponíveis; não aplicar limite superior comum entre plataformas;
4. excluir genes detetados em menos de 1% das células malignas ou menos de 20 células;
5. `log1p(counts/library_size × 10 000)`;
6. selecionar 2 000 genes de dispersão residual dentro de cada coorte, sem outcomes, forçando apenas a interseção final para matching;
7. remover para fitting genes de ciclo celular, stress técnico/heat-shock, mitocondriais, ribossomais e imunoglobulinas; estes programas são pontuados depois como potenciais confundidores;
8. amostrar no máximo 2 000 células por doente para impedir domínio por profundidade; seed 20260906.

As decisões e contagens antes/depois de QC serão registadas por coorte. Não haverá integração batch-corrected antes da descoberta separada.

## Descoberta intra-coorte

- NMF não negativo em expressão log-normalizada não centrada.
- Ranks candidatos 5–15; seeds fixas 11, 23, 47, 71 e 101.
- Escolha de rank por cophenetic/stabilidade entre seeds, erro de reconstrução e ausência de fragmentação extrema; sem outcomes.
- Consenso entre seeds por matching ótimo de correlação Spearman de pesos.
- Um candidato exige estabilidade mediana entre seeds ≥0,60 e pelo menos 20 genes com peso positivo relevante.
- Programas dominados por ciclo, stress técnico, mitocondriais, ribossomais ou imunoglobulinas são rotulados como provável artefacto/confundidor, mas não apagados silenciosamente.

Para GSE246613, os nove programas já congelados são a representação imutável dessa coorte; não se repete NMF. Todos os nove entram no matching com igual estatuto. P8 não recebe tratamento especial.

## Correspondência inter-coorte

Para cada par de candidatos, calcular em genes comuns:

- correlação Spearman dos pesos completos;
- Jaccard dos 50 genes principais;
- overlap ponderado dos 100 genes principais;
- correlação de scores em pseudobulk/células apenas como diagnóstico, nunca entre células não correspondentes.

Correspondência primária exige simultaneamente Spearman ≥0,30 e Jaccard top-50 ≥0,15, com ambos acima do percentil 95 de 10 000 pares de assinaturas permutadas e matched por expressão. Arestas são construídas sem nomes biológicos; famílias são componentes obtidas após matching máximo um-para-um por par de coortes.

## Presença, ausência e indeterminação

Para cada `State × Dataset`:

- **PRESENTE:** candidato correspondente passa os dois limiares, null empírico e estabilidade intra-coorte;
- **AUSENTE:** cobertura ≥80% dos top-50 do estado, ≥1 000 células malignas e nenhum candidato passa, com limite superior bootstrap abaixo de pelo menos um limiar primário;
- **INDETERMINADO:** cobertura/profundidade insuficiente, instabilidade, resultados discordantes ou ausência sem potência suficiente.

Ausência não é inferida apenas porque o matching não é significativo.

## Critério de atlas e confiança

- **Recorrente:** presente em ≥2 coortes de construção.
- **Robusto:** presente em ≥3 coortes de construção, incluindo ≥2 centros/plataformas, estabilidade mediana ≥0,70 e não dominado por confundidor.
- **Validado:** presente na coorte Pal reservada usando assinatura e limiares congelados.
- **Possivelmente universal:** robusto e validado, presente em ≥75% das coortes elegíveis; este rótulo não significa universalidade biológica literal.

Confiança alta requer validação reservada; confiança moderada requer ≥3 coortes de construção; restantes estados são baixa confiança ou específicos de coorte.

## Estabilidade e sensibilidades

- bootstrap de doentes (1 000) e subsampling de células a 25%, 50% e 75%;
- repetição sem os 10 genes de maior peso;
- features 1 500/2 500 e ranks adjacentes como sensibilidades, sem escolher o resultado favorável;
- leave-one-cohort-out no vocabulário final;
- estratificação pós-definição por plataforma, tratamento e baseline/on-treatment;
- uma célula nunca é unidade de validação independente.

## Hierarquia

Será explorada apenas após congelar as famílias de estados, por distância entre assinaturas e co-ocorrência paciente-nível. Uma árvore/rede é descritiva; não implica transição, ancestralidade ou trajetória.

## Relação com P8, CNV e outcomes

- P8 é comparado ao atlas apenas no fim, com os mesmos critérios de matching.
- CNV pode caracterizar estados depois de estabelecidos; não define estados nesta fase.
- pCR, RD, R2/NR e sobrevivência não participam em discovery, rank, matching ou validação de existência.

## Critério de decisão

A hipótese de um atlas comum é suportada apenas se existir pelo menos um estado robusto e validado e se a estrutura global exceder os nulls de permutação. Se nenhum estado validar, concluir: **OS DADOS PÚBLICOS ATUAIS NÃO SUPORTAM UM ATLAS UNIVERSAL DE ESTADOS MALIGNOS.**

Este documento será acompanhado por um SHA-256 antes do primeiro novo download da Fase 3.
