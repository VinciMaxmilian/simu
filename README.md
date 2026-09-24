# simu

Civilizações simuladas pesquisando um universo regido pelo Modelo Padrão.

Quanto tempo leva para surgir uma civilização científica? O que culturas de pesquisa diferentes
descobrem quando a física do laboratório é **calculada** (não sorteada)? Alguma delas acha algo
permitido pelo Modelo Padrão que nós ainda não fizemos?

Design completo em [`docs/PLANO.md`](docs/PLANO.md); dados do Modelo Padrão em
[`docs/pesquisa/modelo_padrao.md`](docs/pesquisa/modelo_padrao.md) e `data/modelo_padrao/*.json`.

## Compilar e rodar

```sh
cmake -S . -B build -G Ninja && cmake --build build

./build/simu emergencia --planetas 20000000                     # camada A: surgimento
./build/simu pesquisa --lab brinquedo --rodadas 32 --inicio 2 --teto 100
./build/simu pesquisa --lab quimica   --rodadas 32 --inicio 5 --teto 500
./build/simu pesquisa --lab materiais --rodadas 32 --inicio 5 --teto 500
./build/simu pesquisa --lab nuclear   --rodadas 32 --inicio 5 --teto 500
./build/simu validar --lab materiais                             # confere contra fatos conhecidos

python3 scripts/graficos.py          # gráficos em resultados/*/grafico.png (requer matplotlib)
python3 scripts/calibrar_nuclear.py  # recalibra o modelo nuclear (requer numpy)
```

Uma rodada completa (as cinco simulações) leva ~2,5 minutos em 4 núcleos.

## Laboratórios

| fase | laboratório | física | validação |
|---|---|---|---|
| 4 | `brinquedo` | paisagem NK com 4 leis escondidas | ótimo global por força bruta (2²⁸) |
| 2 | `quimica` | Hückel em benzenoides; modos zero = magnetismo de carbono | benzeno…[3]triangulena exatos |
| 1 | `materiais` | tight-binding 2D; bandas planas com **gap certificado** | quadrada, Lieb, kagome; estresse sem violações |
| 3 | `nuclear` | gota líquida + camadas calibradas; α, β, fissão | AME2020/NUBASE (qualitativo) |

## Resultados

Cada `resultados/<lab>/resumo.md` traz validação, desempenho por cultura, catálogo de descobertas
e a história da ciência de um universo. Os CSVs completos ficam ao lado.
