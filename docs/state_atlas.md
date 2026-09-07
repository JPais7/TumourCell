# Atlas computacional de estados malignos do TNBC — versão 1

Data: 2026-09-06  
Âmbito: configurações transcricionais malignas; definição independente de outcomes  
Estado: atlas de construção congelado, seguido de validação reservada negativa/inconclusiva

## Resultado executivo

Foram identificadas quatro famílias recorrentes entre coortes de construção. Apenas `State_01` cumpriu a definição pré-especificada de robustez na construção (três coortes e duas plataformas). Nenhuma família foi validada na coorte Pal reservada. Por conseguinte, este documento é um vocabulário **provisório de baixa confiança**, não um atlas universal.

As frequências celulares não são reportadas como proporções nesta versão: os estados foram definidos como programas contínuos de coexpressão e os artefactos congelados não contêm uma calibração comum de uso celular. Uma frequência aparente entre plataformas seria enganadora. A recorrência abaixo refere-se a coortes, não a células.

## Vocabulário provisório

| Estado | Assinatura resumida | Construção | Estabilidade mediana | Pal | Interpretação possível, não mecanística | Limitação principal |
|---|---|---:|---:|---|---|---|
| `State_01` | MGP, VIM, TUBA1A, FOS, CALD1, TACSTD2, JUNB, KRT19, CCND1, SFRP1 | ARTEMIS, Wu, Karaayvaz | 0,865 | indeterminado | eixo epitelial–mesenquimal com resposta imediata/stress | mistura de biologia e dissociação; não observado em GSE246613 e não validado |
| `State_02` | S100A6, HLA-B, HLA-A, RARRES1, KRT19, SAA1, IFI27, TACSTD2, S100A16, CD74 | GSE246613-P8, ARTEMIS, Wu | 0,865 | indeterminado | apresentação antigénica/interferão com componente epitelial | apenas 10x na construção; possível dependência de plataforma; não validado |
| `State_03` | HLA-B, IFI27, HLA-A, IFITM3, LDHB, S100A6, CD74, VIM, BST2, MT2A | ARTEMIS, Wu | 0,859 | indeterminado | variante interferão/HLA | duas coortes; possível fragmentação de `State_02` |
| `State_04` | MUCL1, AZGP1, NEAT1, JUN, KRT19, JUNB, FOS, KLF6, SELENOP, IER2 | ARTEMIS, Wu | 0,859 | indeterminado | diferenciação epitelial/luminal e resposta imediata | duas coortes; pode refletir composição, plataforma ou estado transitório |

As assinaturas completas e os membros de cada família estão em `results/phase3/state_atlas_v1.json` e `results/state_metadata.csv`.

## Hierarquia emergente

Os dados não suportam uma árvore de estados. `State_02` e `State_03` partilham HLA-A/HLA-B, IFI27, CD74, BST2 e IFITM3, sendo compatíveis com variantes ou extremos de um eixo comum de interferão/apresentação antigénica. Separá-los como estados universais distintos é **NÃO SUPORTADO**. `State_01` e `State_04` combinam marcadores de identidade epitelial com genes de resposta imediata; não é possível distinguir estado biológico, dissociação ou mistura dos dois sem validação ortogonal.

## Relação com P8

P8 entrou no matching com o mesmo estatuto dos restantes nove programas congelados de GSE246613. Correspondia a `State_02` na construção através de ARTEMIS e Wu. Contudo, `State_02` não cumpriu o requisito de duas plataformas e não validou no Pal. Logo, P8 não é resgatado como biomarcador nem como estado universal; a conclusão das Fases 1–2 permanece inalterada.

## O que este atlas permite — e não permite

Permite formular assinaturas candidatas e testes de replicação futuros. Não permite inferir transições, clones, resposta, prognóstico, causalidade, dinâmica ou um Digital Twin. A associação Estado → Clone por CNV e a validação espacial ficam condicionadas ao estabelecimento futuro de pelo menos um estado em holdout independente.

**Conclusão:** **OS DADOS PÚBLICOS ATUAIS NÃO SUPORTAM UM ATLAS UNIVERSAL DE ESTADOS MALIGNOS.**
