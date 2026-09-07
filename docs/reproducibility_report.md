# Relatório de reprodutibilidade — Fase 3

Data: 2026-09-06

## Desenho e prevenção de leakage

Os critérios, papéis das coortes, limiares e regras de decisão foram congelados antes dos novos downloads (`results/phase3_protocol_freeze.sha256`). GSE246613, ARTEMIS, Wu/GSE176078 e Karaayvaz/GSE118389 foram atribuídos à construção; Pal/GSE161529 foi reservado. Outcomes não foram usados na seleção de features, rank, descoberta, matching ou validação. O atlas de construção foi congelado antes da abertura analítica do Pal (`results/phase3_construction_freeze.sha256`). Os quatro ficheiros desse manifesto continuam a passar SHA-256.

## Dados efetivamente analisados

| Coorte | Papel | Células malignas usadas | Doentes | Plataforma/matriz | Rank | Estabilidade mediana |
|---|---|---:|---:|---|---:|---:|
| GSE246613 | construção, programas previamente congelados | conforme Fase 1 | conforme Fase 1 | 10x | 9 programas | congelada |
| ARTEMIS | construção | 45 051 de 49 275 | 97 | 10x, counts em H5AD | 5 | 0,853 |
| Wu/GSE176078 | construção | 8 649 amostradas de 10 836 | 8 | 10x, counts em H5AD | 8 | 0,865 |
| Karaayvaz/GSE118389 | construção | 1 169 de 1 534 | 6 | Fluidigm C1, RSEM expected counts | 5 | 0,881 |
| Pal/GSE161529 | holdout | 12 262 amostradas de 36 226 pós-QC | 8 | 10x, counts em H5AD integrado | 9 | 0,825 |

O equilíbrio por doente limitou a 2 000 células/doente. No Pal, um doente tinha apenas 51 células após os filtros; a heterogeneidade de profundidade reduz a potência para declarar ausência.

## Descoberta e matching de construção

A NMF foi executada separadamente em cada coorte, com ranks 5–15, cinco seeds fixas e inicialização aleatória durante 20 épocas. O rank escolhido foi o menor a até 1% da maior estabilidade mediana. Não houve integração batch-corrected antes da descoberta.

O matching exigiu simultaneamente Spearman ≥0,30, Jaccard top-50 ≥0,15 e ambos acima do percentil 95 de 10 000 nulls estratificados por expressão. Das correspondências ótimas, seis arestas passaram: quatro ARTEMIS–Wu, uma Wu–Karaayvaz e uma GSE246613–ARTEMIS. Produziram quatro componentes recorrentes; apenas `State_01` atravessou três coortes e duas plataformas.

## Validação reservada Pal

O Pal foi processado sem modificar o atlas. Para chamar `PRESENTE`, um único candidato Pal tinha de corresponder independentemente a pelo menos dois membros congelados do estado. Resultado:

| Estado | Correspondências Pal que passaram | Estado final |
|---|---:|---|
| `State_01` | 0 | INDETERMINADO |
| `State_02` | 1 (ARTEMIS) | INDETERMINADO |
| `State_03` | 1 (ARTEMIS) | INDETERMINADO |
| `State_04` | 1 (ARTEMIS) | INDETERMINADO |

Houve quatro correspondências significativas ARTEMIS–Pal, incluindo uma assinatura ARTEMIS não pertencente ao atlas, mas nenhuma confirmação convergente por um segundo membro congelado. O padrão é compatível com efeito de coorte/plataforma ou com fragmentação/rotação dos componentes NMF. Não é prova de ausência biológica.

## Robustez e limitações

- Estabilidade intra-coorte: adequada em mediana (0,825–0,881), mas quantis inferiores foram mais fracos, sobretudo Karaayvaz e Pal.
- Estabilidade inter-coorte: limitada; apenas uma família atingiu três coortes/duas plataformas na construção e nenhuma validou.
- Plataforma: `State_02`–`State_04` são suportados apenas por dados 10x; dependência de plataforma é possível e **NÃO IDENTIFICÁVEL** com segurança. `State_01` inclui Fluidigm C1, mas falhou no holdout 10x.
- Tratamento: ARTEMIS, Wu, Karaayvaz e Pal são essencialmente baseline/primário, enquanto GSE246613 tem contexto terapêutico distinto. Só `State_02` inclui GSE246613; não existem réplicas suficientes por tratamento para estimar dependência. **NÃO IDENTIFICÁVEL**.
- Número de células: o cap por doente reduziu domínio de grandes coortes, mas a assimetria 1 169–45 051 células permanece. Não foi estimado um efeito causal do tamanho amostral. **NÃO IDENTIFICÁVEL**.
- Frequência celular comparável: **NÃO IDENTIFICÁVEL**, porque os componentes são contínuos e não existe limiar calibrado comum congelado.
- A auditoria 3b executou 20 bootstraps de doentes, 15 subsamplings por coorte (25%, 50% e 75% por doente), ranks adjacentes e features 1 500/2 500. A remoção dos dez genes de maior peso foi aplicada às arestas do atlas. Estes ensaios são diagnósticos Monte Carlo, não uma nova descoberta nem os 1 000 bootstraps confirmatórios do protocolo completo. `State_01` foi estável em ARTEMIS, moderado em Wu e frágil em Karaayvaz (recuperação mediana ~0,53–0,55); Wu mostrou fragilidade adicional ao bootstrap em `State_02`/`State_04`. As correlações das arestas permaneceram >0,30 após remoção dos genes dominantes, mas isso não restaura a validação no Pal.
- A análise 3b reforça que a ausência de matching não deve ser declarada `AUSENTE`: há instabilidade paciente/plataforma, sobretudo na coorte Fluidigm pequena, e a confiança final permanece baixa.
- Karaayvaz disponibiliza RSEM expected counts fracionários, não UMI counts; a comparação é tecnicamente mais distante, mas também menos diretamente comparável.
- O H5AD Pal provém de uma integração epitelial pública; foram usados apenas os oito identificadores TNBC/BRCA1-TNBC recuperados do estudo original.

## Reprodutibilidade computacional

Os inputs têm SHA-256 registado; scripts, seeds, ranks e limiares estão versionados no repositório. `src/phase3_discover_states.py` e `src/phase3_match_states.py` mantêm os hashes do congelamento de construção. O PDF foi inspecionado e representa apenas a rede congelada de construção; o holdout não altera as arestas.

Resultado global: há recorrência parcial acima dos nulls par-a-par, mas não validação independente suficiente para uma estrutura universal. **NÃO SUPORTADO**.
