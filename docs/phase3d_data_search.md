# Fase 3d — Procura dirigida de dados adicionais

Data da auditoria: 8 de setembro de 2026.

## Objetivo

Procurar dados capazes de resolver as três lacunas deixadas pela Fase 3c:

1. uma segunda validação clínica independente de pCR/resposta;
2. amostragem longitudinal de células malignas TNBC;
3. validação espacial com cobertura génica suficiente.

Foram privilegiados dados processados, públicos, com identificação de doente, momento de colheita e endpoint clínico. Não foram descarregados ficheiros da ordem dos GB.

## Resultado principal

O melhor novo conjunto é o **NeoTRIP, GSE319641**. O ensaio randomizou 280 doentes TNBC de alto risco para carboplatina+nab-paclitaxel, com ou sem atezolizumab. A expressão processada contém 401 biópsias de 251 doentes: 241 basais, 160 em D1C2 e 150 pares completos. A matriz TPM/ComBat pública ocupa cerca de 44 MB.

A cobertura é adequada para projeção: 100/100 genes de cada um dos quatro estados congelados e 1848/2000 genes de P8 estão presentes. Assim, este conjunto permite testar imediatamente mudanças precoces dos estados e de P8 dentro do mesmo doente.

Contudo, o GEO contém apenas doente, momento e lote. O braço terapêutico e o pCR não são públicos; o repositório analítico oficial indica que os metadados clínicos devem ser pedidos a `translational@fondazionemichelangelo.org`. Sem esses campos, o NeoTRIP não pode ainda substituir uma validação clínica de pCR, embora seja já utilizável para uma análise longitudinal não supervisionada.

## Prioridades executáveis

### 1. GSE319641 / NeoTRIP

- Fazer projeção congelada de P8 e State_01–State_04 em todas as 401 amostras.
- Testar deltas basal→D1C2 nos 150 pares, primeiro globalmente e depois por padrões individuais.
- Manter qualquer associação a resposta ou interação com atezolizumab bloqueada até obter os metadados clínicos.
- Dados já descarregados: matriz de expressão de 44 MB e matriz GEO de 14 KB.

### 2. GSE260693

- 67 doentes TNBC; 42 biópsias pré-tratamento e 31 resseções pós-tratamento.
- Inclui anotação pública de tumor residual e identificadores que permitem análise longitudinal.
- É o melhor candidato imediatamente utilizável para replicar associação basal com resposta e estudar persistência/seleção após tratamento.
- Os ficheiros TPM e counts têm menos de 11 MB cada; apenas a matriz de metadados foi descarregada nesta auditoria.

### 3. BioKey / EGAS00001004809

- Estudo single-cell pre/on-treatment com anti-PD1, incluindo TNBC e células malignas.
- O laboratório disponibiliza publicamente contagens por doente com metadados de pareamento, tratamento, expansão imunitária e estado hormonal; os dados brutos permanecem controlados na EGA.
- É a melhor oportunidade identificada para dinâmica maligna single-cell, mas é necessário confirmar quantos pares TNBC contêm células malignas suficientes antes de o promover a validação formal.

### 4. HRA002137

- 62 578 células de sete doentes TNBC HER2-low/negative, amostrados antes, durante e depois de terapêutica neoadjuvante.
- Muito relevante biologicamente, mas de acesso controlado e pequeno para inferência clínica definitiva.

## Coortes clínicas adicionais

- **GSE25066:** grande coorte neoadjuvante com pCR/RD; útil após filtragem rigorosa para TNBC, embora use microarrays antigos.
- **GSE18864:** 24 TNBC antes de cisplatina, com resposta Miller–Payne; serve como teste de sensibilidade específico ao regime, não como validação principal.
- **BrighTNess:** 634 TNBC randomizados, RNA-seq basal e pCR; seria excelente, mas não foi localizada expressão individual aberta.
- **GeparNuevo:** existem resultados públicos para 2549 genes e scripts; os dados clínicos/individuais do ensaio requerem candidatura à GBG.

## Exclusões relevantes

- **GSE266919:** 44 doentes e 78 biópsias, mas os objetos processados públicos contêm apenas células imunitárias (B, CD4, CD8, mieloides e NK). Não valida estados malignos e sobrepõe-se ao contexto do GSE169246.
- **GSE289825:** apenas dois doentes; a série tumoral longitudinal existe só no não respondedor e usa painel imunitário dirigido.
- **GSE296517:** CosMx 1000-plex, essencialmente dois doentes do ensaio; objetos de 1.1–3.3 GB e raw de 175.6 GB. Não foi descarregado.
- **GSE325216:** série espacial de quatro doentes metastáticos, mas só um é TNBC e o painel tem 960 genes.

## Lacuna espacial

Não foi encontrado um conjunto espacial TNBC independente que reúna simultaneamente: vários doentes, amostras clinicamente informativas, expressão whole-transcriptome e acesso público manejável. O Visium HD **EGAD50000002284** tem melhor potencial de cobertura, mas é controlado e inclui um BAM de 60.6 GB. A validação espacial formal deve permanecer assinalada como lacuna, não ser simulada com painéis inadequados.

## Decisão recomendada

Avançar em duas vias:

1. analisar já o GSE260693 para resposta/resíduo e o GSE319641 para dinâmica longitudinal congelada;
2. pedir os metadados clínicos do NeoTRIP e auditar/descarregar as contagens processadas BioKey apenas depois de confirmar o tamanho e a composição TNBC maligna.

Esta combinação acrescenta evidência independente sem inflacionar a confiança: GSE260693 responde à questão clínica disponível agora; NeoTRIP responde à dinâmica precoce; BioKey poderá responder à dinâmica ao nível de célula maligna.
