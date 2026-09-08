# Fase 3d — Validação clínica e dinâmica longitudinal em bulk RNA-seq

Data: 8 de setembro de 2026.

## Perguntas

1. P8 e os quatro estados congelados distinguem resposta patológica numa segunda coorte TNBC independente?
2. Estes programas mudam durante o tratamento no mesmo tumor?

As definições de P8 e State_01–State_04 foram mantidas sem refitting, seleção de genes por outcome ou otimização de limiar.

## GSE260693: validação clínica independente

A expressão pública contém 73 amostras de 54 doentes únicos: 42 pré-tratamento, 31 pós-tratamento e 19 pares completos. Entre as 42 amostras basais, 16 pertencem a doentes sem tumor residual e 26 a doentes com tumor residual. A análise usou log1p(TPM).

A cobertura é completa para os quatro estados (100/100 genes em cada) e quase completa para P8 (1886/2000 genes; 97.47% do peso congelado).

### Resultado basal

O contraste foi orientado como tumor residual `sim` menos `não`. Valores negativos significam, portanto, maior expressão no grupo sem tumor residual.

| Programa | Hedges g, residual−sem residual | IC bootstrap da diferença | p | FDR |
|---|---:|---:|---:|---:|
| State_01 | -0.101 | [-0.121, 0.085] | 0.766 | 0.948 |
| State_02 | -0.158 | [-0.152, 0.087] | 0.386 | 0.643 |
| State_03 | -0.329 | [-0.181, 0.053] | 0.228 | 0.571 |
| State_04 | -0.007 | [-0.122, 0.120] | 0.948 | 0.948 |
| P8 | -0.431 | [-0.282, 0.037] | 0.117 | 0.571 |

P8 e State_03 replicam a direção favorável observada no I-SPY2, mas o GSE260693 isolado não tem precisão suficiente para uma conclusão formal.

## Meta-análise I-SPY2 + GSE260693

Os efeitos foram orientados como resultado favorável (pCR/sem tumor residual) menos doença residual. Foi usado um modelo de efeitos fixos porque existem apenas duas coortes; os efeitos individuais e a heterogeneidade são sempre reportados.

| Programa | g I-SPY2 | g GSE260693 | g combinado [IC95%] | p | FDR | I² |
|---|---:|---:|---:|---:|---:|---:|
| State_01 | 0.176 | 0.101 | 0.168 [-0.032, 0.369] | 0.0993 | 0.1241 | 0% |
| State_02 | 0.260 | 0.158 | 0.249 [0.048, 0.450] | 0.0150 | 0.0281 | 0% |
| State_03 | 0.328 | 0.329 | 0.328 [0.127, 0.529] | 0.00138 | 0.00689 | 0% |
| State_04 | -0.059 | 0.007 | -0.052 [-0.252, 0.148] | 0.609 | 0.609 | 0% |
| P8 | 0.223 | 0.431 | 0.245 [0.044, 0.445] | 0.0169 | 0.0281 | 0% |

O sinal mais convincente é **State_03**: os efeitos das duas coortes são praticamente idênticos. **P8** tem concordância direcional e associação combinada, mas a validação externa isolada continua inconclusiva; deve ser descrito como suporte adicional, não como biomarcador validado.

## Dinâmica no GSE260693

Nos 19 pares pré→pós-tratamento, P8 diminuiu: delta médio -0.307, IC bootstrap [-0.512, -0.114], Wilcoxon p=0.00618, FDR=0.0309. Nenhum estado passou FDR; State_04 mostrou apenas tendência de aumento (delta médio 0.109, FDR=0.136).

As amostras cirúrgicas pós-tratamento requerem tumor disponível e representam sobretudo doença residual. Esta análise informa a biologia das células/tecidos persistentes e não deve ser interpretada como evolução de toda a população tratada.

## Dinâmica precoce no NeoTRIP / GSE319641

A análise incluiu 401 amostras de 251 doentes: 241 basais, 160 em D1C2 e 150 pares completos. A cobertura foi 100/100 genes em cada estado e 1848/2000 genes de P8, correspondendo a 96.76% do peso.

| Programa | Delta médio D1C2−basal | IC bootstrap | Wilcoxon p | FDR |
|---|---:|---:|---:|---:|
| State_01 | -4.464 | [-7.993, -0.982] | 0.0374 | 0.0623 |
| State_02 | 1.338 | [-3.549, 5.970] | 0.239 | 0.299 |
| State_03 | 4.351 | [-0.383, 9.041] | 0.0212 | 0.0529 |
| State_04 | -5.457 | [-9.755, -1.195] | 0.00185 | 0.00924 |
| P8 | 1.199 | [-4.513, 6.937] | 0.303 | 0.303 |

O resultado robusto é a redução precoce de **State_04**. State_03 aumenta e State_01 diminui, mas ficam imediatamente acima do limiar FDR. P8 não apresenta alteração média precoce detectável.

Sem o braço terapêutico e o pCR, não é possível determinar se estas alterações dependem de atezolizumab ou predizem resposta. Esses metadados têm de ser solicitados aos autores.

## Conclusão

- **State_03:** melhor candidato clínico atual; replicação direcional exata e associação combinada robusta.
- **P8:** apoio clínico adicional e redução em doença residual após tratamento, mas ainda não atinge validação externa autónoma.
- **State_04:** programa mais claramente modulado no início do NeoTRIP, independentemente do outcome ainda desconhecido.
- **Dados em falta:** braço e pCR do NeoTRIP continuam a ser o passo com maior valor marginal.

## Ficheiros reproduzíveis

- `src/analyze_phase3d_bulk.py`
- `results/phase3d/gse260693_frozen_scores.csv`
- `results/phase3d/gse260693_validation.json`
- `results/phase3d/neotrip_frozen_scores.csv`
- `results/phase3d/neotrip_longitudinal.json`
- `results/phase3d/ispy2_gse260693_meta_analysis.json`
