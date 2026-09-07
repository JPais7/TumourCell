# Robustez e influência

## População principal

O limiar congelado de 50 células malignas por momento produz 12 doentes: 6 NR, 5 R2 e 1 R1. Consequentemente, R1 é apenas descritivo. A inferência entre grupos concentra-se em R2 versus NR.

## Resultado candidato P8

P8, apresentação antigénica/interferão, foi o único programa com contraste R2–NR positivo nas duas etapas:

| Limiar | Etapa | Efeito R2–NR (SD) | IC bootstrap 95% |
|---:|---|---:|---:|
| 25 | pembrolizumab | 0,268 | −0,139 a 0,690 |
| 25 | adição de radioterapia | 0,268 | −0,527 a 0,996 |
| 50 | pembrolizumab | 0,343 | −0,078 a 0,781 |
| 50 | adição de radioterapia | 0,226 | −0,626 a 1,006 |
| 100 | pembrolizumab | 0,270 | −0,120 a 0,619 |
| 100 | adição de radioterapia | 0,342 | −0,072 a 0,736 |

Nenhum intervalo exclui zero. O sinal é consistente em direção entre ranks 8–10, cinco seeds, 1 500/2 000/2 500 genes, três normalizações e limiares 25/50/100. Depois de ajuste por P2/P3, o efeito primário da segunda etapa é 0,302 SD. Esta consistência não compensa a baixa precisão.

## Influência

No contraste P8 da segunda etapa com limiar 50, a direção permanece em 10 de 11 análises leave-one-out. A remoção de Patient45 inverte o efeito para −0,029 SD; Patient63 e Patient19 também alteram materialmente a magnitude. No limiar 100 a direção permanece estável, mas restam apenas 3 R2 e 4 NR.

## Alterações longitudinais não específicas da resposta

P3 diminui após a adição de radioterapia tanto em R2 como NR. P2 e P4 diminuem em NR. Estes resultados são exploratórios, não corrigidos para multiplicidade e parcialmente correlacionados com expressão total; são mais compatíveis com efeitos gerais de tratamento/amostragem do que com um marcador robusto de resposta.

## Explicações alternativas

- Profundidade associa-se moderada/fortemente a vários deltas na segunda etapa.
- Batch está confounded com biopsia/momento e não pode ser separado plenamente.
- P4 sugere possível contaminação immune-like na anotação maligna.
- Os scores ponderados correlacionam com expressão total selecionada.
- A composição de estados malignos pode mudar sem alteração dentro das células.
- R1 tem n=1 no conjunto primário.

## Conclusão

A direção de P8 é metodologicamente estável, mas estatisticamente inconclusiva e influenciável. Nenhum programa satisfaz um critério forte de associação diferencial com resposta.
