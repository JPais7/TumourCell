# Descoberta cega de programas malignos

Data: 2026-09-06  
Estado: concluída e congelada antes da análise de resposta

## Separação metodológica

O script de descoberta não carrega `response_group`, sobrevivência, doença residual ou qualquer endpoint. A seleção de genes, a amostragem, os ranks e os seeds dependem apenas de expressão, anotação maligna, doente e momento. A resposta foi aberta apenas pelo processo posterior `analyze_frozen_programs.py`, depois de verificado o hash da definição congelada.

## Input e seleção de genes

- População: 50 185 células com `subtype_new == cancer cells`.
- Matriz: `X`, expressão não negativa normalizada/log-transformada; os counts de `raw/X` ficaram reservados para pseudobulk.
- Deteção exigida: 1%–90% das células malignas.
- Excluídos: genes `MT-`, `RPL*` e `RPS*`.
- Mantidos: genes de ciclo celular e stress, para identificação e ajuste explícitos.
- Ranking: variância logarítmica padronizada dentro de 20 bins de expressão média.
- Input primário: 2 000 genes.
- Amostragem balanceada: até 100 células por doente×momento, seed técnico 246613.
- Escala: divisão pelo desvio-padrão do gene nas células malignas, sem centrar.

## Comparação de ranks

Foram testados ranks 4–10 e seeds 11, 29, 47, 71 e 101.

| Rank | Erro médio | Estabilidade média | Estabilidade mínima | Redundância máxima média |
|---:|---:|---:|---:|---:|
| 4 | 3536,56 | 0,954 | 0,770 | 0,529 |
| 5 | 3521,32 | 0,949 | 0,871 | 0,617 |
| 6 | 3506,78 | 0,936 | 0,893 | 0,624 |
| 7 | 3494,23 | 0,948 | 0,870 | 0,669 |
| 8 | 3483,15 | 0,961 | 0,885 | 0,664 |
| 9 | 3472,55 | 0,960 | 0,913 | 0,660 |
| 10 | 3463,80 | 0,946 | 0,830 | 0,684 |

O rank 9 foi escolhido por apresentar a melhor estabilidade mínima, estabilidade média quase máxima e convergência nos cinco seeds. Um modelo rank 8 atingiu o limite de 500 iterações. O seed 11 foi o medoide por similaridade média aos outros seeds do rank 9.

## Artefactos

- `results/phase1/blind_discovery/gene_selection.json`
- `results/phase1/blind_discovery/rank_seed_metrics.json`
- `results/phase1/blind_discovery/candidate_program_weights.npz`
- `figures/phase1/nmf_rank_stability.png`

Nenhuma comparação R1/R2/NR influenciou esta escolha.
