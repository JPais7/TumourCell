# Validação externa

## Conclusão

**VALIDAÇÃO EXTERNA NÃO DISPONÍVEL** para a pergunta longitudinal maligna completa sob pembrolizumab seguido de pembrolizumab+radioterapia.

## Melhor teste parcial disponível

Yan/ARTEMIS 2026 é o conjunto independente mais informativo para testar se o score baseline congelado de P8 se associa a pCR versus doença residual. Não pode validar a alteração longitudinal nem distinguir seleção de plasticidade, porque o componente molecular relevante é pré-tratamento.

Antes de abrir ARTEMIS ficam congelados:

- os 2 000 genes e pesos `program_definitions_v1`;
- score contínuo sem threshold;
- hipótese direcional exploratória: maior P8 baseline associa-se a resposta favorável;
- unidade de inferência: doente;
- sem retuning, nova seleção de genes ou novo rank.

## Critérios

- **Replicação parcial:** associação de P8 baseline na direção prevista, com efeito não dominado por um doente e estabilidade a profundidade/qualidade.
- **Falha:** efeito consistentemente oposto com incerteza suficientemente estreita.
- **Inconclusivo:** intervalo amplo, poucos doentes elegíveis, transferência incompleta de genes ou dependência de poucos doentes.

GSE169246 não serve para esta validação porque contém células CD45+ e não o compartimento maligno. Kim/SRP114962 seria a extensão mais informativa para seleção/plasticidade se a matriz processada e o dicionário celular se tornarem auditáveis.
