# Relatório científico da Fase 1

## Resultado principal

A análise identificou nove programas reprodutíveis de expressão maligna. O sinal candidato mais coerente com resposta tardia é P8, apresentação antigénica/interferão, mas o efeito R2–NR após adição de radioterapia é pequeno e impreciso: 0,226 SD, IC bootstrap 95% −0,626 a 1,006. Nenhum programa diferencia robustamente R2 de NR.

## Programas e estabilidade

O rank 9, seed 11, foi congelado sem acesso à resposta. A estabilidade média entre seeds foi 0,960 e a mínima 0,913. Os programas abrangem housekeeping/citoesqueleto, stress, ciclo celular, expressão immune-like/CD45, estados androgénico e apócrino, basal/mesenquimal, apresentação antigénica/interferão e luminal/secretor.

## Alterações longitudinais

P3 diminuiu depois da radioterapia em R2 e NR, sugerindo redução proliferativa não específica da classe de resposta. P2 e P4 diminuíram em NR. P8 manteve direção positiva R2–NR em todas as sensibilidades, mas todos os intervalos incluíram zero e Patient45 consegue inverter a direção no leave-one-out primário.

## Seleção e plasticidade

Todos os programas são mensuráveis no baseline, mas isto não demonstra populações preexistentes discretas. Não existem CNV/clones no objeto e as células não estão ligadas entre tempos.

**INDETERMINADO — os dados atuais não permitem distinguir seleção de reprogramação transcricional.**

## Robustez

- Direção P8 preservada em ranks 8–10 e seeds alternativos.
- Direção preservada com 1 500 e 2 500 genes.
- Direção preservada nos limiares de 25, 50 e 100 células.
- Direção preservada em log-CPM, log-TP10k e sqrt-CPM.
- Efeito P8 após ajuste por stress/ciclo: 0,302 SD.
- Precisão insuficiente em todas as variantes.

## Conclusões por nível epistemológico

### OBSERVADO

- Nove eixos de covariação reproduzíveis na expressão maligna.
- Alterações longitudinais de score dentro de doentes.
- Redução de P3 após a segunda etapa em R2 e NR.

### INFERIDO

- P8 representa variação relacionada com apresentação antigénica/interferão.
- Parte substancial da variação pode refletir expressão total, composição e qualidade de amostragem.

### SUPORTADO

- O modelo NMF é estável entre seeds.
- Não existe evidência robusta de que um programa diferencie R2 de NR nesta amostra.

### INCERTO

- P8 pode acompanhar a resposta tardia, mas a amostra é demasiado pequena para uma conclusão firme.
- A relevância do padrão observado no único R1 elegível.

### NÃO IDENTIFICÁVEL

- Seleção versus plasticidade.
- Persistência ou expansão clonal.
- Causalidade da radioterapia.

## Validação e próximo experimento mínimo

Não existe validação externa longitudinal adequada. ARTEMIS pode testar apenas a associação baseline de P8 com resposta, usando a assinatura já congelada. O próximo experimento computacional mínimo para abordar seleção/plasticidade é inferir CNV por doente e momento com referência normal adequada, validar a qualidade dessas chamadas e testar alterações de programa dentro de perfis copy-number persistentes. Se a resolução CNV falhar, a conclusão deve permanecer não identificável.

## Decisão

**PARAR a escalada de complexidade.** Não construir digital twin, classificador clínico ou modelo de transição. O único avanço justificado é validação externa parcial de P8 e uma auditoria técnica de viabilidade de CNV inferido.
