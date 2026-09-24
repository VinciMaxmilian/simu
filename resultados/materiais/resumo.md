# Materiais: redes 2D com bandas planas (tight-binding)

**Objetivo de pesquisa:** achar a banda mais plana possível e, de preferência, isolada das outras por um gap de energia; score = −log10(largura + 1e−6) + 3·min(gap certificado, 1).

32 universos × 4 civilizações × 300 anos de ciência; 7198925 experimentos (1678273 distintos) em 38.3 s.

## Validação contra o mundo real

| rede | esperado (teoria conhecida) | modelo |
|---|---|---|
| quadrada | uma banda, E ∈ [−4, 4], largura 8 | banda mais plana: E=0.000, largura=8.00e+00, gap certificado ≥ 0.000 (aparente 0.000); faixa total [-4.000, 4.000] → `dispersiva` |
| Lieb (CuO₂ dos cupratos) | banda plana em E=0 tocando as outras (gap 0) → não isolada | banda mais plana: E=0.000, largura=3.99e-16, gap certificado ≥ 0.000 (aparente 0.000); faixa total [-2.828, 2.828] → `banda plana nao isolada [desbalanco de subredes], 3 sitios` |
| kagome | banda plana em E=+2|t| tocando em Γ (gap 0) → não isolada | banda mais plana: E=2.000, largura=4.00e-15, gap certificado ≥ 0.000 (aparente 0.000); faixa total [-4.000, 2.000] → `banda plana nao isolada [interferencia], 3 sitios` |
| falso isolado da v1, analisado na grade 12×12 onde a v1 errava | banda plana tocando (as bandas se encontram fora da grade 12×12) | banda mais plana: E=0.000, largura=2.03e-16, gap certificado ≥ 0.000 (aparente 1.000); faixa total [-2.828, 2.828] → `banda plana nao isolada [desbalanco de subredes], 3 sitios` |
| isolada achada pelas civilizações | banda plana isolada (gap ≈ 1.236 na grade 24×24) | banda mais plana: E=-0.000, largura=5.87e-16, gap certificado ≥ 0.989 (aparente 1.236); faixa total [-3.236, 3.236] → `banda plana isolada [interferencia], 3 sitios` |

Teste de estresse do certificado: 300 redes aleatórias com banda (quase) plana e gap certificado, 600000 pontos k contínuos sorteados, 0 violações (menor folga observada: 0.001).

## Desempenho por cultura de pesquisa

| cultura | explora | recombina | contato | melhor score (média ± dp) | classes descobertas | descobertas inéditas no universo | ano médio em que atingiu o recorde global | catástrofes |
|---|---|---|---|---|---|---|---|---|
| Empiristas | 0.50 | 0.10 | 0.15 | 9 ± 1.1e-11 | 8.0 | 81 | 58 (32/32 rodadas) | 1.2 |
| Incrementalistas | 0.03 | 0.05 | 0.15 | 9 ± 1.1e-11 | 8.0 | 59 | 57 (32/32 rodadas) | 1.1 |
| Evolucionistas | 0.10 | 0.45 | 0.15 | 9 ± 1.1e-11 | 8.0 | 53 | 56 (32/32 rodadas) | 1.6 |
| Isolados | 0.15 | 0.20 | 0.00 | 9 ± 1.2e-11 | 7.4 | 68 | 78 (32/32 rodadas) | 1.2 |

## Catálogo de descobertas (11 classes; ordenadas pelo melhor score)

| classe | melhor score | 1º ano (qualquer universo) | universos que acharam |
|---|---|---|---|
| banda plana isolada [desbalanco de subredes], 3 sitios | 9 | 2 | 32/32 |
| banda plana isolada [interferencia], 3 sitios | 9 | 2 | 32/32 |
| banda plana isolada [interferencia], 4 sitios | 9 | 1 | 32/32 |
| banda plana nao isolada [interferencia], 4 sitios | 6.149 | 7 | 32/32 |
| banda plana nao isolada [desbalanco de subredes], 3 sitios | 6.143 | 1 | 32/32 |
| banda plana nao isolada [interferencia], 3 sitios | 6.141 | 0 | 32/32 |
| banda plana nao isolada [desbalanco de subredes], 4 sitios | 6 | 0 | 30/32 |
| banda plana nao isolada [interferencia], 2 sitios | 6 | 264 | 1/32 |
| quase plana isolada, 4 sitios | 2.352 | 169 | 2/32 |
| quase plana nao isolada, 4 sitios | 1.647 | 168 | 4/32 |
| dispersiva | 0.4235 | 0 | 32/32 |

## Melhores exemplares

**banda plana isolada [desbalanco de subredes], 3 sitios** (score 9.0)

```
3 sítios/célula; banda 2 de 3 em E=-0.000|t|, largura=3.16e-16 (grade 24×24), gap certificado ≥ 0.997 (aparente na grade: 1.247)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 0→1(0,0)+1 0→1(1,0)-1 1→0(1,0)-1 0→1(0,1)-1 0→2(0,1)-1 1→0(0,1)-1 0→1(1,1)+1 0→2(1,1)+1 1→0(1,1)-1 0→1(1,-1)-1
```

**banda plana isolada [interferencia], 3 sitios** (score 9.0)

```
3 sítios/célula; banda 2 de 3 em E=0.000|t|, largura=6.76e-16 (grade 24×24), gap certificado ≥ 0.986 (aparente na grade: 1.232)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 0→2(0,0)+1 0→2(1,0)-1 2→0(1,0)-1 2→2(1,0)+1 2→1(0,1)-1 0→2(1,1)+1 2→0(1,1)-1 0→2(1,-1)-1 1→2(1,-1)-1
```

**banda plana isolada [interferencia], 4 sitios** (score 9.0)

```
4 sítios/célula; banda 2 de 4 em E=-1.000|t|, largura=3.22e-15 (grade 24×24), gap certificado ≥ 1.131 (aparente na grade: 1.414)
    energias locais: 0 0 -1 0
    saltos (i→j, célula, t): 0→1(0,0)+1 2→1(1,0)-1 3→2(1,0)+1 0→2(0,1)-1 2→1(1,1)+1 1→2(1,-1)-1 3→0(1,-1)+1 3→1(1,-1)+1
```

**banda plana nao isolada [interferencia], 4 sitios** (score 6.1)

```
4 sítios/célula; banda 2 de 4 em E=-1.000|t|, largura=4.22e-15 (grade 24×24), gap certificado ≥ 0.071 (aparente na grade: 0.089)
    energias locais: 0 0 0 -1
    saltos (i→j, célula, t): 0→1(0,0)-1 1→3(0,0)-1 0→2(1,0)-1 1→0(1,0)+1 1→2(1,0)+1 3→1(1,0)-1 1→0(0,1)-1 1→2(0,1)+1 3→1(1,1)+1 1→0(1,-1)-1 2→1(1,-1)-1 3→1(1,-1)-1
```

**banda plana nao isolada [desbalanco de subredes], 3 sitios** (score 6.1)

```
3 sítios/célula; banda 2 de 3 em E=-0.000|t|, largura=3.94e-16 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.436)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 2→0(1,0)-1 0→2(0,1)-1 1→2(0,1)+1 2→0(0,1)+1 2→1(1,1)-1 0→2(1,-1)-1 2→0(1,-1)-1
```

**banda plana nao isolada [interferencia], 3 sitios** (score 6.1)

```
3 sítios/célula; banda 2 de 3 em E=-0.000|t|, largura=6.61e-16 (grade 24×24), gap certificado ≥ 0.041 (aparente na grade: 0.082)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 0→1(1,0)-1 0→2(1,0)+1 1→0(1,0)+1 0→1(0,1)+1 1→0(0,1)-1 2→0(1,1)+1 0→0(1,-1)+1 0→1(1,-1)-1 0→2(1,-1)-1 1→0(1,-1)+1
```

**banda plana nao isolada [desbalanco de subredes], 4 sitios** (score 6.0)

```
4 sítios/célula; banda 2 de 4 em E=-0.000|t|, largura=3.90e-16 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.000)
    energias locais: 0 0 0 0
    saltos (i→j, célula, t): 0→1(0,0)+1 0→1(1,0)-1 1→0(1,0)-1 3→0(1,0)+1 0→1(0,1)-1 0→2(0,1)-1 0→1(1,1)+1 0→2(1,1)+1 1→0(1,1)-1 0→3(1,-1)-1
```

**banda plana nao isolada [interferencia], 2 sitios** (score 6.0)

```
2 sítios/célula; banda 2 de 2 em E=2.000|t|, largura=1.33e-15 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.000)
    energias locais: 0 0
    saltos (i→j, célula, t): 0→1(0,0)-1 0→1(1,0)-1 0→0(0,1)-1 1→0(0,1)-1 0→1(1,1)-1 1→1(1,1)-1
```

## História da ciência no universo 0 (primeiras descobertas)

- ano 0: **Empiristas** descobrem `dispersiva`
- ano 5: **Isolados** descobrem `banda plana isolada [interferencia], 4 sitios`
- ano 13: **Empiristas** descobrem `banda plana isolada [interferencia], 3 sitios`
- ano 17: **Evolucionistas** descobrem `banda plana isolada [desbalanco de subredes], 3 sitios`
- ano 19: **Isolados** descobrem `banda plana nao isolada [interferencia], 4 sitios`
- ano 22: **Incrementalistas** descobrem `banda plana nao isolada [interferencia], 3 sitios`
- ano 22: **Incrementalistas** descobrem `banda plana nao isolada [desbalanco de subredes], 3 sitios`
- ano 203: **Incrementalistas** descobrem `banda plana nao isolada [desbalanco de subredes], 4 sitios`
