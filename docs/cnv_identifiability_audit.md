# Auditoria de identificabilidade de CNV e clonalidade — Fase 2

Data da auditoria: 2026-09-06  
Âmbito: GSE246613 e Yan/ARTEMIS; distinção entre seleção e plasticidade

## Conclusão executiva

**SELECTION VS PLASTICITY: NÃO IDENTIFICÁVEL COM OS DADOS ATUAIS.**

Os dados de scRNA-seq permitem procurar padrões de expressão ordenados genomicamente e, em ARTEMIS, os autores publicaram uma classificação consensual de células aneuploides. Isso é adequado como apoio à identificação de células malignas. Não equivale a medir CNV em DNA, definir clones biologicamente defensáveis ou demonstrar que um clone no baseline é o mesmo clone numa biopsia residual. ARTEMIS é, além disso, um atlas pré-tratamento: não contém no objeto público analisado a observação longitudinal necessária para correspondência clonal antes/depois do tratamento.

## Evidência técnica auditada

### GSE246613

- O objeto disponível é scRNA-seq; a matriz bruta contém counts de expressão.
- Não foram encontrados perfis CNV, segmentos, genótipos somáticos, scDNA-seq, WES/WGS pareado ou identificadores clonais validados no objeto utilizado na Fase 1.
- Existem três momentos de amostragem, mas biopsias separadas não permitem rastrear diretamente a mesma célula.
- É tecnicamente possível inferir sinal CNV-like a partir de RNA, mas a referência normal, cobertura génica, dropout, estados de expressão e composição de cada biopsia condicionam fortemente o resultado.

### ARTEMIS

- O objeto público CELLxGENE é scRNA-seq de tecido não tratado, com uma anotação final de células anormais/aneuploides.
- O código dos autores corre CopyKAT por amostra misturando células normais mamárias externas (epiteliais, imunes e estromais) com as células observadas.
- A decisão publicada combina: classificação hierárquica CopyKAT, magnitude/correlação CNA tipo Tirosh e clustering Leiden com normais como controlo negativo. Uma quarta estratégia baseada em TCGA foi conservada no código, mas não usada no manuscrito.
- Parâmetros registados no código: hg38, janela mínima de 25 genes, `KS.cut=0.2`; o resultado é uma inferência baseada em expressão, não uma medição de DNA.
- Os próprios autores avisam que sobre-expressão de genes vizinhos — explicitamente genes HLA no cromossoma 6 — pode produzir eventos CNA focais espúrios, e que o baseline da referência altera a inferência.

## Cadeia de inferência

| Nível | Afirmação | Dados necessários | Dados disponíveis / método possível | Resolução, validação e confundidores | Veredicto |
| --- | --- | --- | --- | --- | --- |
| 1 | Sinal de expressão compatível com CNV | Counts de RNA, ordem genómica e referência | Disponível em ambas as coortes; smoothing por posição/CopyKAT ou método equivalente | Resolução baixa; confundido por programas coordenados, profundidade, dropout, ciclo, stress, tipo celular e batch | **Possível, exploratório** |
| 2 | Segmentos CNV-like inferidos | Sinal nível 1 consistente em janelas e referência apropriada | Possível; ARTEMIS usa janelas de pelo menos 25 genes | Segmentos não são chamadas de DNA; focalidade e limites são instáveis. HLA/IFN é um confundidor direto para P8 | **Possível, baixa confiança** |
| 3 | Perfis CNV por célula/doente | Cobertura suficiente e chamadas estáveis por célula | ARTEMIS disponibiliza a classificação final aneuploide; perfis RNA-CNV podem ser inferidos | Dropout torna perfis celulares ruidosos; consenso por doente melhora estabilidade, mas reduz resolução subclonal | **Parcialmente possível** |
| 4 | Subpopulações com perfis CNV distintos | Perfis estáveis, distância/threshold congelados e separação robusta | Clustering exploratório seria executável | Clusters dependem do método, referência e threshold; falta validação ortogonal e estabilidade demonstrada entre biopsias | **Não estabelecido** |
| 5 | Clones biologicamente defensáveis | Alterações de DNA partilhadas, idealmente SNV/CNV ortogonais e filogenia | Não disponíveis nos objetos scRNA-seq auditados | Um cluster RNA-CNV não identifica necessariamente uma linhagem clonal | **Não identificável** |
| 6 | Correspondência de clones entre momentos/biopsias | Marcadores clonais de DNA partilhados e amostras longitudinais pareadas | GSE246613 é longitudinal mas sem marcadores clonais validados; o objeto ARTEMIS é baseline | CNV convergente, ruído e amostragem espacial impedem identidade clonal segura | **Não identificável** |

Alcançar os níveis 1–3 não autoriza as conclusões dos níveis 4–6.

## Relação entre P8 e CNV

P8 é dominado por apresentação antigénica/interferão, incluindo vários genes HLA. Um aumento coordenado desses genes pode parecer uma alteração focal do cromossoma 6 num método RNA-CNV. Assim:

1. a anotação aneuploide publicada pode ser usada como regra independente de inclusão de células malignas;
2. não se deve usar um segmento chr6/HLA inferido da mesma expressão para “validar” P8 — seria circular;
3. associação entre P8 e score RNA-CNV não demonstra que P8 é clonal;
4. diferenças de P8 entre biopsias não demonstram, isoladamente, plasticidade.

## Hipóteses concorrentes

| Hipótese | O que os dados podem testar nesta fase | O que permanece em falta |
| --- | --- | --- |
| H1 — associação biológica real | Validação externa paciente-nível do P8 congelado em ARTEMIS | Mecanismo e causalidade |
| H2 — seleção | Enriquecimento fenotípico pode ser descrito | Clone preexistente rastreável e sua expansão |
| H3 — plasticidade | Mudança longitudinal de score pode ser descrita em GSE246613 | Identidade da mesma linhagem/célula e exclusão de seleção/composição |
| H4 — composição | Sensibilidades ao número de células, pseudobulk e anotação maligna | Amostragem espacial completa da lesão |
| H5 — efeito técnico | Sensibilidades de normalização, profundidade e batch | Eliminação completa de todos os confundidores técnicos |
| H6 — especificidade da coorte | Teste externo sem retuning | Mais coortes independentes e endpoint plenamente homólogo |

## Validação necessária para ultrapassar o limite

Para atingir níveis 5–6 seriam necessários, no mínimo, dados pareados baseline/residual por doente com scDNA-seq, WES/WGS ou painel de DNA suficientemente informativo; alterações somáticas partilhadas para definir clones; controlo de pureza e copy number; e concordância entre DNA e RNA. Idealmente, códigos de barras lineage-tracing seriam necessários num sistema experimental para distinguir diretamente persistência clonal de mudança de estado.

## Decisão

- **CNV-like a partir de RNA:** identificável apenas de forma exploratória (níveis 1–3).
- **Clonalidade biologicamente defensável:** **NÃO IDENTIFICÁVEL**.
- **Clone baseline → mesmo clone residual:** **NÃO IDENTIFICÁVEL**.
- **Seleção versus plasticidade:** **NÃO IDENTIFICÁVEL COM OS DADOS ATUAIS**.

Não será construída uma narrativa de expansão clonal a partir de clusters RNA-CNV, nem uma narrativa de plasticidade a partir de alterações de score.

## Proveniência

- Repositório ARTEMIS auditado: `navinlabcode/tnbc-chemo`, commit `fbe1dd3b1db05dd5e9f02485a45fdd438d6af9fb`.
- Script dos autores: `analysis/scripts/copykat_mix.R`.
- Nota metodológica dos autores: `analysis/identifying_aneuploid_cells.md`.
- Objeto público: `e94bd3cc-6271-424a-baac-12f8eb320a0e`, coleção `ceef2841-5333-46ac-92ef-ccbe0c20fe55`.
