# Mundo-brinquedo (paisagem NK com leis escondidas)

**Objetivo de pesquisa:** maximizar o desempenho de uma 'tecnologia' de 28 componentes; score = fração do ótimo global.

32 universos × 4 civilizações × 300 anos de ciência; 1866310 experimentos (897703 distintos) em 0.6 s.

## Validação contra o mundo real

Ótimo global (força bruta sobre 2^28 configurações): config 0001001001011011101101011010 — 100.00% do ótimo
Leis satisfeitas pelo ótimo global: Lei-A+Lei-C+Lei-D
Lei-A: x8=0 x19=1 x2=0 x13=0 x11=1 (bônus 0.12)
Lei-B: x19=0 x8=0 x7=1 x16=1 x18=1 (bônus 0.12)
Lei-C: x25=0 x10=0 x17=0 x27=0 x8=0 (bônus 0.12)
Lei-D: x19=1 x26=1 x9=1 x24=1 x16=1 (bônus 0.12)

## Desempenho por cultura de pesquisa

| cultura | explora | recombina | contato | melhor score (média ± dp) | classes descobertas | descobertas inéditas no universo | ano médio em que atingiu o recorde global | catástrofes |
|---|---|---|---|---|---|---|---|---|
| Empiristas | 0.50 | 0.10 | 0.15 | 0.9902 ± 0.011 | 10.0 | 141 | 171 (13/32 rodadas) | 0.9 |
| Incrementalistas | 0.03 | 0.05 | 0.15 | 0.9903 ± 0.011 | 10.0 | 78 | 172 (13/32 rodadas) | 1.4 |
| Evolucionistas | 0.10 | 0.45 | 0.15 | 0.9907 ± 0.011 | 10.0 | 54 | 178 (14/32 rodadas) | 1.2 |
| Isolados | 0.15 | 0.20 | 0.00 | 0.9892 ± 0.012 | 9.9 | 47 | 235 (12/32 rodadas) | 0.8 |

## Catálogo de descobertas (10 classes; ordenadas pelo melhor score)

| classe | melhor score | 1º ano (qualquer universo) | universos que acharam |
|---|---|---|---|
| Lei-A+Lei-C+Lei-D | 1 | 3 | 32/32 |
| Lei-B+Lei-C | 0.9878 | 0 | 32/32 |
| Lei-C+Lei-D | 0.9081 | 6 | 32/32 |
| Lei-A+Lei-C | 0.8962 | 2 | 32/32 |
| Lei-A+Lei-D | 0.8919 | 2 | 32/32 |
| Lei-B | 0.8612 | 0 | 32/32 |
| Lei-C | 0.8582 | 0 | 32/32 |
| Lei-A | 0.8395 | 0 | 32/32 |
| Lei-D | 0.7983 | 0 | 32/32 |
| nenhuma lei | 0.7109 | 0 | 32/32 |

## Melhores exemplares

**Lei-A+Lei-C+Lei-D** (score 1.0)

```
config 0001001001011011101101011010 — 100.00% do ótimo
```

**Lei-B+Lei-C** (score 1.0)

```
config 0001101100010001101000010000 — 98.78% do ótimo
```

**Lei-C+Lei-D** (score 0.9)

```
config 0011110001001111101110001010 — 90.81% do ótimo
```

**Lei-A+Lei-C** (score 0.9)

```
config 0001101100010001101100010000 — 89.62% do ótimo
```

**Lei-A+Lei-D** (score 0.9)

```
config 1001110101010011111110011111 — 89.19% do ótimo
```

**Lei-B** (score 0.9)

```
config 1001110100011001101000000100 — 86.12% do ótimo
```

**Lei-C** (score 0.9)

```
config 0001101101000011100000010000 — 85.82% do ótimo
```

**Lei-A** (score 0.8)

```
config 0000101001010011110101010001 — 83.95% do ótimo
```

## História da ciência no universo 0 (primeiras descobertas)

- ano 0: **Empiristas** descobrem `nenhuma lei`
- ano 0: **Isolados** descobrem `Lei-B`
- ano 2: **Isolados** descobrem `Lei-C`
- ano 3: **Empiristas** descobrem `Lei-A`
- ano 26: **Isolados** descobrem `Lei-D`
- ano 32: **Empiristas** descobrem `Lei-A+Lei-C`
- ano 35: **Evolucionistas** descobrem `Lei-B+Lei-C`
- ano 63: **Isolados** descobrem `Lei-A+Lei-C+Lei-D`
- ano 82: **Incrementalistas** descobrem `Lei-A+Lei-D`
- ano 109: **Isolados** descobrem `Lei-C+Lei-D`
