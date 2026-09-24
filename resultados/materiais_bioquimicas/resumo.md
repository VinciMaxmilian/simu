# Materiais: redes 2D com bandas planas (tight-binding)

**Objetivo de pesquisa:** achar a banda mais plana possível e, de preferência, isolada das outras por um gap de energia; score = −log10(largura + 1e−6) + 3·min(gap certificado, 1).

32 universos × 6 civilizações × 300 anos de ciência; 7203332 experimentos (2147255 distintos) em 123.1 s.

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
| Carbono/agua | 0.15 | 0.20 | 0.00 | 9 ± 1.2e-11 | 7.7 | 52 | 85 (32/32 rodadas) | 1.3 |
| Carbono/amonia | 0.15 | 0.20 | 0.00 | 9 ± 1.1e-11 | 7.5 | 34 | 99 (32/32 rodadas) | 1.2 |
| Carbono/acido sulfurico | 0.15 | 0.20 | 0.00 | 9 ± 1.1e-11 | 7.5 | 82 | 65 (32/32 rodadas) | 1.4 |
| Carbono/metano | 0.15 | 0.20 | 0.00 | 9 ± 5.6e-11 | 6.6 | 15 | 116 (32/32 rodadas) | 1.6 |
| Silicio/agua | 0.15 | 0.20 | 0.00 | 9 ± 1.1e-11 | 7.6 | 38 | 81 (32/32 rodadas) | 1.0 |
| Silicio/acido sulfurico | 0.15 | 0.20 | 0.00 | 9 ± 1.2e-11 | 7.5 | 55 | 70 (32/32 rodadas) | 1.3 |

## Catálogo de descobertas (11 classes; ordenadas pelo melhor score)

| classe | melhor score | 1º ano (qualquer universo) | universos que acharam |
|---|---|---|---|
| banda plana isolada [desbalanco de subredes], 3 sitios | 9 | 1 | 32/32 |
| banda plana isolada [interferencia], 4 sitios | 9 | 1 | 32/32 |
| banda plana isolada [interferencia], 3 sitios | 9 | 2 | 32/32 |
| banda plana nao isolada [interferencia], 4 sitios | 6.146 | 2 | 32/32 |
| banda plana nao isolada [interferencia], 3 sitios | 6.143 | 0 | 32/32 |
| banda plana nao isolada [desbalanco de subredes], 3 sitios | 6 | 0 | 32/32 |
| banda plana nao isolada [desbalanco de subredes], 4 sitios | 6 | 0 | 30/32 |
| banda plana nao isolada [interferencia], 2 sitios | 6 | 217 | 1/32 |
| quase plana isolada, 4 sitios | 4.302 | 64 | 8/32 |
| quase plana nao isolada, 4 sitios | 2.587 | 24 | 13/32 |
| dispersiva | 0.4781 | 0 | 32/32 |

## Melhores exemplares

**banda plana isolada [desbalanco de subredes], 3 sitios** (score 9.0)

```
3 sítios/célula; banda 2 de 3 em E=0.000|t|, largura=4.05e-16 (grade 24×24), gap certificado ≥ 1.386 (aparente na grade: 1.732)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 0→1(0,1)-1 1→0(0,1)-1 0→1(1,1)-1 1→0(1,1)+1 1→2(1,1)+1 2→1(1,-1)-1
```

**banda plana isolada [interferencia], 4 sitios** (score 9.0)

```
4 sítios/célula; banda 2 de 4 em E=-0.000|t|, largura=3.05e-16 (grade 24×24), gap certificado ≥ 1.131 (aparente na grade: 1.414)
    energias locais: 0 1 1 0
    saltos (i→j, célula, t): 0→3(0,0)-1 0→3(1,0)-1 3→0(1,0)+1 3→1(1,0)+1 0→3(0,1)-1 2→1(0,1)-1 3→0(0,1)+1 0→3(1,-1)+1 3→2(1,-1)+1
```

**banda plana isolada [interferencia], 3 sitios** (score 9.0)

```
3 sítios/célula; banda 2 de 3 em E=0.000|t|, largura=6.62e-16 (grade 24×24), gap certificado ≥ 1.042 (aparente na grade: 1.303)
    energias locais: 0 0 1
    saltos (i→j, célula, t): 0→2(0,1)-1 2→0(0,1)-1 1→2(1,1)+1 0→2(1,-1)-1 2→0(1,-1)+1 2→1(1,-1)-1
```

**banda plana nao isolada [interferencia], 4 sitios** (score 6.1)

```
4 sítios/célula; banda 3 de 4 em E=0.000|t|, largura=1.08e-15 (grade 24×24), gap certificado ≥ 0.052 (aparente na grade: 0.103)
    energias locais: 0 0 -1 0
    saltos (i→j, célula, t): 1→2(0,0)+1 1→0(1,0)+1 1→1(1,0)+1 1→1(0,1)-1 1→2(0,1)-1 1→3(0,1)-1 3→1(0,1)+1 1→0(1,-1)-1 1→3(1,-1)+1
```

**banda plana nao isolada [interferencia], 3 sitios** (score 6.1)

```
3 sítios/célula; banda 2 de 3 em E=0.000|t|, largura=9.06e-16 (grade 24×24), gap certificado ≥ 0.094 (aparente na grade: 0.118)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 0→1(0,0)-1 0→0(1,0)+1 0→1(0,1)-1 0→2(0,1)-1 1→0(0,1)-1 2→0(1,1)-1 0→0(1,-1)-1 0→1(1,-1)-1 1→0(1,-1)-1 2→0(1,-1)-1
```

**banda plana nao isolada [desbalanco de subredes], 3 sitios** (score 6.0)

```
3 sítios/célula; banda 2 de 3 em E=0.000|t|, largura=3.80e-16 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.000)
    energias locais: 0 0 0
    saltos (i→j, célula, t): 1→2(1,0)-1 1→0(0,1)-1 1→2(0,1)+1 1→0(1,1)-1 1→0(1,-1)-1
```

**banda plana nao isolada [desbalanco de subredes], 4 sitios** (score 6.0)

```
4 sítios/célula; banda 3 de 4 em E=0.000|t|, largura=4.59e-16 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.000)
    energias locais: 0 0 0 0
    saltos (i→j, célula, t): 2→0(1,0)+1 2→1(1,0)-1 2→0(0,1)-1 2→3(0,1)-1 0→2(1,1)-1 2→1(1,1)-1 2→3(1,1)-1 2→0(1,-1)+1
```

**banda plana nao isolada [interferencia], 2 sitios** (score 6.0)

```
2 sítios/célula; banda 1 de 2 em E=-2.000|t|, largura=1.11e-15 (grade 24×24), gap certificado ≥ 0.000 (aparente na grade: 0.000)
    energias locais: 0 0
    saltos (i→j, célula, t): 0→1(0,0)-1 0→1(1,0)-1 0→1(0,1)+1 1→1(0,1)-1 0→0(1,-1)-1 0→1(1,-1)+1
```

## História da ciência no universo 0 (primeiras descobertas)

- ano 0: **Carbono/agua** descobrem `dispersiva`
- ano 2: **Silicio/agua** descobrem `banda plana nao isolada [interferencia], 3 sitios`
- ano 5: **Silicio/agua** descobrem `banda plana nao isolada [desbalanco de subredes], 3 sitios`
- ano 6: **Carbono/metano** descobrem `banda plana nao isolada [interferencia], 4 sitios`
- ano 6: **Silicio/agua** descobrem `banda plana isolada [desbalanco de subredes], 3 sitios`
- ano 12: **Carbono/metano** descobrem `banda plana isolada [interferencia], 4 sitios`
- ano 18: **Silicio/acido sulfurico** descobrem `banda plana isolada [interferencia], 3 sitios`
- ano 105: **Carbono/acido sulfurico** descobrem `quase plana isolada, 4 sitios`
- ano 113: **Carbono/acido sulfurico** descobrem `quase plana nao isolada, 4 sitios`
- ano 151: **Carbono/metano** descobrem `banda plana nao isolada [desbalanco de subredes], 4 sitios`
