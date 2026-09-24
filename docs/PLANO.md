# simu — design

## Pergunta

Num universo regido pelo Modelo Padrão, quanto tempo leva para surgir uma civilização consciente,
quanto tempo ela leva para fazer ciência, o que ela pesquisa — e ela acharia algo, **permitido pelo
Modelo Padrão**, que nós ainda não fizemos?

## Por que não simular átomo por átomo

| alvo | custo aproximado |
|---|---|
| recorde atual de dinâmica molecular | ~10¹² átomos por nanossegundos |
| 1 s de um cérebro humano (~10²⁶ átomos, passo 1 fs) | ~10⁴¹ passos-átomo ≈ 10¹⁸ anos em exaescala |
| 4 bilhões de anos de Terra (~10⁵⁰ átomos) | ~10⁸² passos-átomo |

Por isso o `simu` é **multiescala**: cada camada usa a física efetiva da sua escala.

## Camadas

**A — Surgimento** (`src/core/emergence.*`). Monte Carlo de planetas: formação ao longo da história
cósmica (Madau–Dickinson, ΛCDM Planck 2018), massa estelar, janela habitável e "passos difíceis"
(abiogênese → eucariontes → multicelularidade → inteligência → ciência).

**B — Civilizações** (`src/core/research.*`). Cada civilização tem uma cultura de pesquisa
(explorar × refinar × recombinar × trocar conhecimento), capacidade experimental crescente,
catástrofes e idades das trevas que apagam parte do conhecimento.

**C — Laboratórios virtuais** (`src/labs/*`). A física real que as civilizações investigam. Para
que uma descoberta delas valha no mundo real, o laboratório tem que **calcular** as consequências
das leis (não sortear um resultado que nós escrevemos). Todo laboratório se valida contra fatos
conhecidos antes de rodar.

| fase | laboratório | física | fidelidade |
|---|---|---|---|
| 4 | brinquedo | paisagem NK com 4 leis escondidas (ótimo global conhecido) | exata (é artificial) |
| 2 | química | Hückel (elétrons π) em benzenoides; modos zero = magnetismo de carbono | boa para essa família |
| 1 | materiais | tight-binding 2D; bandas planas | exata dentro do modelo |
| 3 | nuclear | gota líquida + camadas calibradas; α, β, fissão | qualitativa |

## Próximos passos

- Química: trocar Hückel por Hubbard/DFT-tight-binding (xTB) para confirmar os candidatos;
  permitir heteroátomos (N, B) e anéis de 5/7 membros.
- Materiais: redes com simetria triangular/hexagonal, orbitais múltiplos, topologia (número de Chern)
  das bandas planas isoladas.
- Nuclear: modelo de massas com deformação (FRDM/Duflo–Zuker) e dados completos do AME2020.
- Ligar A → B: civilizações começam a ciência no instante sorteado pela camada A, com capacidade
  experimental dependente do planeta.
- Civilizações com agentes de IA propondo experimentos (em vez de busca estocástica).
