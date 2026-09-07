# Especificação bloqueável da análise de Fase 1

Estado: protocolo v0.2 congelado antes da análise de expressão

## Pergunta primária

No GSE246613, que programas malignos mudam entre baseline, pós-pembrolizumab e pós-pembrolizumab+radioterapia, e que alterações distinguem resposta precoce (R1), resposta tardia após radioterapia (R2) e ausência de resposta (NR)?

Esta é uma pergunta de redistribuição longitudinal de programas, não uma alegação de plasticidade dentro de clones. Sem lineage tracing ou CNA longitudinal resolvido, seleção e reprogramação permanecem mecanismos compatíveis, mas não identificáveis separadamente.

## Unidade estatística

O doente é a unidade independente. Células e biopsias técnicas são medições aninhadas. O doente 03 do GSE246613 conta uma vez, apesar de ter duas unidades tumorais (`Patient03T1` e `Patient03T2`).

## Coortes e funções bloqueadas

- **Descoberta longitudinal primária:** GSE246613, restrita a células anotadas como `cancer cells`.
- **Extensão clone-aware futura:** Kim 2018, apenas se for recuperada uma matriz processada com doente, tempo e evidência de malignidade/CNA.
- **Validação de associação pré-tratamento:** Yan/ARTEMIS 2026, completamente vedado à definição do programa.
- **Validação externa pré-tratamento:** Yan/ARTEMIS 2026 permanece vedado durante a descoberta e só poderá receber assinaturas congeladas.
- **Validação mecanística imunitária:** GSE169246, paclitaxel ou paclitaxel+atezolizumab; não utilizar para descoberta ou validação de programas malignos porque o perfil é CD45+.
- **Controlo de subtipo:** GSE205472 apenas exploratório, devido aos cinco pares e metadados incompletos.

## Representação de estado

Utilizar programas contínuos de expressão, não clusters exclusivos. Os programas serão definidos por NMF/cNMF nas células malignas do GSE246613, exigindo estabilidade de genes e scores entre seeds e reamostragens de doentes. Ciclo celular, interferão, hipoxia, apoptose, resposta a dano e stress dissociativo serão programas explícitos de controlo. O número de programas será escolhido por estabilidade/reconstrução antes de testar associação com resposta.

## População maligna no GSE246613

O objeto não imune contém 171 746 células, das quais 50 185 estão anotadas como `cancer cells`. Existem 34 doentes e 35 unidades tumorais; resposta por doente: 9 R1, 14 R2 e 11 NR.

Elegibilidade longitudinal primária proposta: pelo menos 50 células malignas em cada um dos três momentos. Isto deixa 12 doentes. Sensibilidades obrigatórias: limiares de 1, 10, 25 e 100 células, que deixam respetivamente 31, 21, 14 e 8 doentes.

O limiar primário fica congelado em 50 células malignas por momento. Os resultados têm de ser apresentados juntamente com todos os limiares de sensibilidade, sem escolher o limiar que produza maior significância.

## Estimandos

1. Alteração intra-doente baseline → pós-pembrolizumab, estratificada por R1/R2/NR.
2. Alteração intra-doente pós-pembrolizumab → pós-pembrolizumab+radioterapia, estratificada por R1/R2/NR.
3. Contraste de interação entre etapa terapêutica e trajetória de resposta, com prioridade para R2 versus NR na segunda etapa.
4. Associação entre o burden pré-tratamento das assinaturas congeladas e pCR/doença residual no ARTEMIS.

## Modelo mínimo

- Pseudobulk por doente/unidade tumoral e tempo; Patient03 será agregado ao nível do doente numa análise e tratado por unidade tumoral numa sensibilidade.
- Contrastes intra-doente com intervalos de confiança por bootstrap de doentes.
- Modelo hierárquico apenas se o número de doentes elegíveis o suportar; a análise principal não depende de valores-p assintóticos com 12 doentes.
- Síntese entre estudos por direção e efeito padronizado; não integrar células numa única matriz como teste primário.

## Falsificadores obrigatórios

- Gene sets aleatórios emparelhados por expressão.
- Downsampling do baseline para a profundidade pós-tratamento e vice-versa.
- Remoção/regressão separada de ciclo, stress, hipoxia e interferão.
- Chamadas alternativas de células malignas e exclusão de células com baixa confiança de anotação.
- Exclusão de cada doente, um de cada vez.
- Teste externo no ARTEMIS sem retuning e teste de especificidade no HR+ quando possível.

## Critério de sucesso de Nível 1

Programa estável em reamostragem de doentes, contraste longitudinal coerente com a etapa terapêutica e trajetória de resposta, e robustez aos principais confundidores. Isto não basta para progressão ao Nível 2.

## Critério de sucesso de Nível 2

Efeito com direção pré-especificada num estudo independente, acompanhado por incerteza útil e sem depender de um único doente. A associação no ARTEMIS tem de ser testada com definições totalmente bloqueadas.

## Critérios de paragem

- Menos de oito doentes são elegíveis depois do critério mínimo de células malignas nos três momentos.
- O programa desaparece após ajuste de stress/ciclo ou muda de sinal entre definições razoáveis.
- A validação externa falha na direção pré-especificada.

## Linguagem permitida

“Associado”, “consistente com”, “preditivo no conjunto externo” e “o modelo prevê”. Não utilizar “causa”, “transição clonal observada”, “plasticidade demonstrada” ou “resistência induzida” sem lineage tracing ou perturbação apropriada.
