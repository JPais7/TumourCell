# Decisão da Fase 3

Data: 2026-09-06  
Decisão: hipótese principal rejeitada no nível de evidência pré-especificado

## Veredicto

**OS DADOS PÚBLICOS ATUAIS NÃO SUPORTAM UM ATLAS UNIVERSAL DE ESTADOS MALIGNOS.**

O critério congelado exigia pelo menos um estado robusto na construção, validado no holdout, e estrutura acima dos nulls. `State_01` foi robusto na construção, mas ficou indeterminado no Pal. Nenhum dos quatro estados recorrentes foi validado. Não se modificou o atlas depois deste resultado.

## Respostas aos critérios de sucesso

1. **Quantos estados malignos reproduzíveis existem?** Quatro famílias recorrentes exploratórias na construção; uma robusta na construção; zero validadas.
2. **Quais aparecem consistentemente entre estudos?** `State_01` aparece em ARTEMIS, Wu e Karaayvaz. `State_02` aparece em GSE246613/P8, ARTEMIS e Wu, mas apenas numa plataforma. `State_03` e `State_04` aparecem em ARTEMIS e Wu.
3. **Quais dependem do tratamento?** **NÃO IDENTIFICÁVEL.** O desenho não tem replicação balanceada por tratamento.
4. **Quais dependem da plataforma?** `State_02`–`State_04` têm suporte restrito a 10x e podem ser dependentes; causalidade de plataforma é **NÃO IDENTIFICÁVEL**. `State_01` atravessa 10x e Fluidigm, mas não validou.
5. **Quais parecem universais?** Nenhum.
6. **Quais são provavelmente artefactos?** Nenhum pode ser declarado artefacto puro. Genes de resposta imediata/stress em `State_01` e `State_04` e a concentração ARTEMIS–Pal sugerem confundimento técnico possível. A distinção biologia versus artefacto é **NÃO IDENTIFICÁVEL** nesta versão.

## Consequências científicas

- P8 mantém o sinal externo fraco/inconclusivo das fases anteriores; a correspondência com `State_02` não o transforma em biomarcador nem em estado universal.
- Seleção versus plasticidade continua não identificável.
- Estado → Clone por CNV não deve avançar como validação central enquanto não existir estado validado.
- Spatial e microambiente podem ser usados futuramente como testes ortogonais, sem causalidade presumida.
- Não há fundamento para transições, ODE, Markov, Neural ODE, perturbações virtuais ou Digital Twin.

## Auditoria de robustez 3b

A auditoria confirmou estabilidade elevada dos programas ARTEMIS, estabilidade intermédia em Wu e fragilidade marcada do membro Karaayvaz de `State_01` (recuperação Spearman mediana 0,53–0,55). A remoção dos dez genes dominantes preservou as correlações das arestas acima de 0,30, sugerindo sinal distribuído, mas não corrigiu a falha de validação no Pal. As sensibilidades de 1 500/2 500 features e ranks adjacentes não transformaram nenhum estado em validado. Portanto, a limitação principal é a generalização entre doentes/coortes, não simplesmente o número total de células.

## Próximo passo defensável

Uma Fase 3c só é justificável com novas coortes TNBC independentes, counts originais, anotação maligna e pelo menos três doentes, aplicando o atlas congelado sem retuning. Idealmente devem incluir ≥20–30 doentes, uma plataforma adicional e dados emparelhados scRNA/snRNA ou spatial. Até essa validação, os quatro estados permanecem candidatos de baixa confiança.
