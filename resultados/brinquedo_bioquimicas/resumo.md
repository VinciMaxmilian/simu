# Mundo-brinquedo (paisagem NK com leis escondidas)

**Objetivo de pesquisa:** maximizar o desempenho de uma 'tecnologia' de 28 componentes; score = fração do ótimo global.

32 universos × 6 civilizações × 300 anos de ciência; 1817853 experimentos (986388 distintos) em 1.0 s.

## Validação contra o mundo real

Ótimo global (força bruta sobre 2^28 configurações): config 0001001001011011101101011010 — 100.00% do ótimo
Leis satisfeitas pelo ótimo global: Lei-A+Lei-C+Lei-D
Lei-A: x8=0 x19=1 x2=0 x13=0 x11=1 (bônus 0.12)
Lei-B: x19=0 x8=0 x7=1 x16=1 x18=1 (bônus 0.12)
Lei-C: x25=0 x10=0 x17=0 x27=0 x8=0 (bônus 0.12)
Lei-D: x19=1 x26=1 x9=1 x24=1 x16=1 (bônus 0.12)

## Desempenho por cultura de pesquisa

Civilizações de bioquímicas diferentes (cenário Terra, catálise adaptada), cada uma no seu sistema estelar (sem contato). Velocidade = pensamento relativo ao humano; ferramentas = fração do teto instrumental (sem fogo/metalurgia: 30%).

| civilização | velocidade | ferramentas |
|---|---|---|
| Carbono/agua | 1.07 | 100% |
| Carbono/amonia | 0.74 | 30% |
| Carbono/acido sulfurico | 1.37 | 30% |
| Carbono/metano | 0.48 | 30% |
| Silicio/agua | 1.11 | 100% |
| Silicio/acido sulfurico | 1.51 | 30% |

| cultura | explora | recombina | contato | melhor score (média ± dp) | classes descobertas | descobertas inéditas no universo | ano médio em que atingiu o recorde global | catástrofes |
|---|---|---|---|---|---|---|---|---|
| Carbono/agua | 0.15 | 0.20 | 0.00 | 0.9881 ± 0.013 | 10.0 | 90 | 190 (11/32 rodadas) | 0.9 |
| Carbono/amonia | 0.15 | 0.20 | 0.00 | 0.9837 ± 0.017 | 9.6 | 36 | 186 (12/32 rodadas) | 1.3 |
| Carbono/acido sulfurico | 0.15 | 0.20 | 0.00 | 0.987 ± 0.015 | 9.9 | 68 | 163 (13/32 rodadas) | 1.3 |
| Carbono/metano | 0.15 | 0.20 | 0.00 | 0.9723 ± 0.02 | 9.3 | 16 | 250 (1/32 rodadas) | 1.2 |
| Silicio/agua | 0.15 | 0.20 | 0.00 | 0.9886 ± 0.013 | 10.0 | 40 | 201 (13/32 rodadas) | 1.1 |
| Silicio/acido sulfurico | 0.15 | 0.20 | 0.00 | 0.9903 ± 0.011 | 9.9 | 70 | 178 (14/32 rodadas) | 1.2 |

## Catálogo de descobertas (10 classes; ordenadas pelo melhor score)

| classe | melhor score | 1º ano (qualquer universo) | universos que acharam |
|---|---|---|---|
| Lei-A+Lei-C+Lei-D | 1 | 10 | 32/32 |
| Lei-B+Lei-C | 0.9878 | 0 | 32/32 |
| Lei-C+Lei-D | 0.9468 | 1 | 32/32 |
| Lei-A+Lei-C | 0.8979 | 2 | 32/32 |
| Lei-A+Lei-D | 0.8953 | 1 | 32/32 |
| Lei-C | 0.8355 | 0 | 32/32 |
| Lei-B | 0.8344 | 0 | 32/32 |
| Lei-D | 0.8292 | 0 | 32/32 |
| Lei-A | 0.8184 | 0 | 32/32 |
| nenhuma lei | 0.6907 | 0 | 32/32 |

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
config 0010110001000111101110001010 — 94.68% do ótimo
```

**Lei-A+Lei-C** (score 0.9)

```
config 0001010100011011001101101000 — 89.79% do ótimo
```

**Lei-A+Lei-D** (score 0.9)

```
config 1000111101010000111100001110 — 89.53% do ótimo
```

**Lei-C** (score 0.8)

```
config 0001101101001011001000010010 — 83.55% do ótimo
```

**Lei-B** (score 0.8)

```
config 0101010100100010101000001101 — 83.44% do ótimo
```

**Lei-D** (score 0.8)

```
config 1011110101010011111110001111 — 82.92% do ótimo
```

## História da ciência no universo 0 (primeiras descobertas)

- ano 0: **Carbono/agua** descobrem `nenhuma lei`
- ano 0: **Silicio/agua** descobrem `Lei-D`
- ano 0: **Silicio/acido sulfurico** descobrem `Lei-B`
- ano 3: **Carbono/agua** descobrem `Lei-C`
- ano 4: **Carbono/agua** descobrem `Lei-C+Lei-D`
- ano 20: **Carbono/agua** descobrem `Lei-A`
- ano 35: **Silicio/agua** descobrem `Lei-A+Lei-D`
- ano 45: **Carbono/acido sulfurico** descobrem `Lei-B+Lei-C`
- ano 48: **Carbono/acido sulfurico** descobrem `Lei-A+Lei-C`
- ano 64: **Carbono/acido sulfurico** descobrem `Lei-A+Lei-C+Lei-D`
