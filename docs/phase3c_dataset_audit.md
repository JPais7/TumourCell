# Fase 3c — Auditoria de coortes adicionais e plano de validação

**Data de corte:** 2026-09-08  
**Âmbito:** validação externa dos estados malignos TNBC já congelados. Esta auditoria não altera a definição dos estados, os genes, os limiares nem os resultados das Fases 1–3b.

## Decisão executiva

1. **HTAN/HTAPP Metastatic Breast Cancer é a prioridade 1.** É a melhor nova coorte para testar generalização para doença metastática, diferentes locais anatómicos e scRNA/snRNA. Deve ser usada como _holdout_ estrito, restringindo a análise principal a biopsias HR−/HER2− e mantendo doente como unidade estatística.
2. **SCP3851 é apenas exploratório para estes estados.** Tem grande número de células malignas e lesões primárias/metastáticas, mas usa um painel Xenium dirigido e as anotações celulares foram transferidas com GSE176078. A auditoria do painel encontrou apenas 22–29/100 genes por estado e 27,2% do peso de P8, abaixo do necessário para validação formal.
3. **I-SPY2 é prioridade 3 para validação clínica ortogonal.** O recurso público principal é bulk de pré-tratamento (microarray/RPPA/imagem), não single-cell. Pode testar associações entre assinaturas congeladas, pCR e braço terapêutico; não demonstra a existência de um estado celular.
4. **GSE176078 e GSE246613 permanecem construção/dinâmica interna.** Já contribuíram para o atlas e não podem ser apresentados como novas validações independentes.
5. **GSE169246 é exclusivamente imune.** A deposição contém células CD45+ e é útil para ligar estados malignos a mecanismos imunes, mas não para validar diretamente estados tumorais.
6. **GSE252315 é excluído como coorte clínica.** É bulk RNA-seq de uma linha HMLE ao longo de TGF-β; serve apenas para triangulação mecanística EMT/plasticidade.

## Auditoria por conjunto pedido

### HTAN/HTAPP Metastatic Breast Cancer — avançar

- 67 biopsias de 60 doentes com cancro da mama metastático; 30 scRNA-seq e 37 snRNA-seq.
- Composição clínica publicada: 44 biopsias HR+/HER2−, 3 HR+/HER2+, 3 HR−/HER2+ e 16 HR−/HER2−; a contagem é de biopsias, não necessariamente de doentes independentes.
- Inclui múltiplos locais metastáticos e até 15 biopsias com modalidades espaciais emparelhadas (Slide-seq, MERFISH, ExSeq, CODEX e H&E).
- O portal disponibiliza contagens, H5AD, tipo celular, compartimento, QC e campo `cnv_pass_mal`. A coleção CELLxGENE auditada nesta data contém 32 objetos, 878 532 células e cerca de 8,4 GB.
- **Risco principal:** forte estrutura por doente, local, subtipo e tecnologia; scRNA e snRNA não devem ser misturados sem estratificação.

**Teste predefinido:** aplicar os programas congelados sem refitting; selecionar HR−/HER2−; confirmar malignidade por anotação e CNV; produzir prevalência e score por biopsia; meta-analisar ao nível do doente; reportar scRNA e snRNA separadamente e em conjunto; repetir por local metastático; excluir doentes/coortes sobrepostos, se existirem.

**Auditoria local do agregado scRNA fresco concluída:** 157 531 células e 26 843 genes, 30 doentes/amostras e 85 033 células malignas. A seleção estrita `receptors_biopsy == ER-/PR-/HER2-` contém 32 205 células de 5 doentes, incluindo 16 724 células malignas. Contagens malignas por doente: 388, 1 131, 2 073, 2 372 e 10 760; todos passam o limiar mínimo. Os quatro estados têm 100/100 genes de topo presentes; P8 tem 1 838/2 000 genes (91,9%). Locais das células malignas TNBC: fígado 10 760, axila 4 445, mama 1 131 e pescoço 388. Este resultado aprova formalmente o agregado fresco para validação. O agregado snRNA congelado deve ser auditado como replicação tecnológica separada.

### GSE176078 — manter como construção

- 26 tumores primários: 11 ER+, 5 HER2+ e 10 TNBC.
- No ficheiro local: 100 064 células no total, 42 512 células TNBC e 10 836 células epiteliais malignas TNBC antes do balanceamento, provenientes de 8 doentes elegíveis.
- Tem excelente utilidade de referência, mas foi usado na construção do atlas; qualquer resultado novo é análise de sensibilidade, não validação independente.

### GSE252315 — excluir da validação clínica

- 21 amostras bulk RNA-seq de HMLE: controlo e TGF-β nos dias 1, 4, 7, 10, 15 e 20, com três réplicas.
- Não contém biopsias de doentes nem scRNA-seq humano. O artigo reutiliza outras coortes single-cell; essas coortes têm de ser avaliadas pelos seus próprios acessos.
- Utilidade legítima: testar se programas congelados seguem a perturbação TGF-β/EMT numa série temporal in vitro, explicitamente rotulada como evidência mecanística.

### GSE169246 — validação imune apenas

- 22 doentes com TNBC avançado, paclitaxel versus paclitaxel + atezolizumab, com amostras pré/pós/progressão.
- Auditoria local: 489 490 células CD45+, 78 amostras RNA (48 sangue, 30 tumor); não existe compartimento maligno adequado.
- Utilidade: relacionar respostas/trajectórias imunes com programas malignos medidos noutras coortes, sem alegar validação direta do atlas tumoral.

### GSE246613 — dinâmica longitudinal já usada

- TNBC antes de tratamento, após pembrolizumab e após pembrolizumab + radioterapia.
- Ficheiro local não imune: 171 746 células, das quais 50 185 anotadas como cancerígenas.
- É forte para trajectórias longitudinais, mas já participou na descoberta/construção. A introdução de radioterapia no terceiro tempo impede atribuir a mudança apenas a pembrolizumab; RNA/CNV inferido também não distingue definitivamente seleção clonal de plasticidade.

### I-SPY2 — validação clínica bulk

- O recurso I-SPY2-990 (GSE194040) inclui expressão pré-tratamento de aproximadamente 19 mil genes em 987 doentes de alto risco, dez braços, subtipo, MammaPrint e pCR; plataforma Agilent, com dados RPPA associados (GSE196093/GSE196096).
- Existem séries específicas de braços, por exemplo GSE173839 e GSE180962, e MRI longitudinal no TCIA.
- **Uso recomendado:** calcular scores pseudobulk usando apenas os genes congelados e testar associação com pCR por subtipo/braço, com modelos que controlem variáveis clínicas e lote. A análise deve ser confirmatória e não usada para redefinir o programa.
- **Limitação:** diferenças entre microarray bulk e scRNA, composição tumoral/TME e ausência de resolução celular.

## Outros dados encontrados

| Prioridade | Conjunto | Valor | Limitação/decisão |
|---|---|---|---|
| 2 | SCP3851, atlas espacial metastático | 126 amostras, 44 doentes, 520 850 células Xenium; 171 766 cancerígenas | Painel 5K; anotação transferida de GSE176078. Validação espacial, não descoberta independente. |
| 4 | GSE208532, derrames pleurais malignos | 10 derrames de 7 doentes metastáticos | Só 2 doentes TNBC; análise exploratória de nicho, abaixo do mínimo para validação. |
| 4 | GSE75688 | Cinco doentes TNBC, tecnologia full-length antiga; inclui uma metástase ganglionar | Cerca de 205 células TNBC: teste de stress de plataforma, não _holdout_ decisivo. |
| acesso controlado | EGAS00001004809, Bassez | Pré/on anti-PD-1, grande atlas single-cell multi-subtipo | Promissor para farmacodinâmica; requer acesso e auditoria do número de TNBC/malignas emparelhadas. |
| acesso controlado | HRA002137 | 62 578 células, 7 TNBC, série durante imunoterapia/NAC | Muito relevante, mas acesso controlado e braços pequenos/heterogéneos. |
| exploratória | GSE264205 | Série temporal e espacial de tumor primário/metástases | Um único doente em estadio IV; estudo de caso, não validação populacional. |
| excluir | GSE138536 | 1 902 células Smart-seq2 de 8 doentes | Só 2 tumores basal-like; seleção de subpopulações epiteliais e desenho de desenvolvimento, insuficiente como TNBC externo. |
| excluir | GSE210616 | TNBC espacial/Visium | Pequena coorte espacial e spots mistos; útil para contexto espacial, não para estados single-cell independentes. |
| mecanística | GSE123837/GSE147326 | Metástase em PDX ou células de efusão cultivadas/perturbadas | Modelos/cultura; triangulação funcional apenas. |

## Critérios de aceitação para uma validação externa

- Estado e pesos génicos congelados antes de abrir o desfecho.
- Pelo menos 3 doentes TNBC e, idealmente, pelo menos 500 células malignas no total e 50 por doente.
- Doente, nunca célula, como unidade de inferência; intervalos por _bootstrap_ de doentes.
- Confirmação de malignidade independente da assinatura testada, preferencialmente anotação original + CNV.
- Separação explícita por tecnologia, tecido/local, tratamento e tempo.
- Resultados negativos e falhas de cobertura génica reportados; sem otimização de limiares na coorte de teste.
- Sobreposição de doentes/amostras com conjuntos de construção verificada antes da análise.

## Ordem de trabalho recomendada

1. Descarregar apenas o H5AD agregado sc/snRNA do HTAN/HTAPP e fazer auditoria de `obs`, genes, doentes, biopsias, subtipo, local, tecnologia, tempo e CNV.
2. Congelar um manifesto HTAN TNBC e uma matriz de testes antes de calcular qualquer score.
3. Executar validação primária por biopsia/doente, com análises separadas scRNA/snRNA e por local. **Concluído para projeções pseudobulk:** 5 doentes/16 724 células malignas scRNA e 11 doentes/68 489 células malignas snRNA.
4. Auditar SCP3851. **Concluído:** cobertura insuficiente (Estado 01 29%, Estado 02 22%, Estado 03 26%, Estado 04 28%; P8 27,2% do peso); não fazer alegação de validação formal.
5. Projetar os programas para I-SPY2 como pseudobulk e testar associação clínica sem retuning. **Concluído:** 362 TNBC (142 pCR, 220 doença residual). Estado 03 associou-se positivamente a pCR (diferença 0,1125; IC bootstrap 0,0433–0,1823; FDR BH 0,0225; coeficiente ajustado por braço/plataforma 0,1018, p=0,0087). Estado 02 teve sinal mais fraco (ajustado p=0,0286; FDR não ajustado 0,0667). P8 foi direcional mas inconclusivo após ajuste (p=0,0596). Estados 01 e 04 não foram confirmados.
6. Usar GSE252315, GSE169246 e coortes PDX/cultura apenas para triangulação mecanística.

## Fontes primárias

- HTAN/HTAPP MBC: DOI 10.1038/s41591-024-03215-z; SCP2702; CELLxGENE collection `a96133de-e951-4e2d-ace6-59db8b3bfb1d`.
- GSE176078: GEO GSE176078; DOI 10.1038/s41588-021-00911-1.
- GSE252315: GEO GSE252315; DOI 10.1038/s44321-024-00050-0.
- GSE169246: GEO GSE169246; DOI 10.1016/j.ccell.2021.09.010.
- GSE246613: GEO GSE246613.
- I-SPY2: GEO GSE194040, GSE196093 e GSE196096.
- SCP3851: Broad Single Cell Portal SCP3851.
- Outros GEO: GSE208532, GSE75688, GSE138536, GSE264205, GSE123837 e GSE147326.
