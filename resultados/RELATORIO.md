# Relatório do simu — 4ª rodada: verificação do C₄₉H₂₁

25 set 2026 · branch `claude/determined-tesla-nkx3fz` · também na `main`

## Resumo

1. **O C₄₉H₂₁ é o único benzenoide plano de até 14 anéis que combina as duas origens de magnetismo**
   (desbalanço de subredes e frustração topológica). A enumeração exaustiva de todos os 15,8 milhões
   de benzenoides de até 14 anéis, validada contra a OEIS, encontra exatamente um, com a mesma
   InChIKey da molécula que as civilizações acharam. Nenhum de até 13 anéis tem essa propriedade.
2. **A busca converge para ele:** foi o resultado final em 265 de 384 civilizações, 86 de 96
   universos, com 3 sementes independentes.
3. **O acoplamento magnético é robusto:** o C₄₉H₂₁ tem estado fundamental de spin 1/2 e um estado de
   spin 3/2 a **~21 meV (faixa 18–22 meV)**. A razão em relação ao cálice de Clar (23 meV medidos)
   fica entre 0,79 e 0,97 em todos os testes: U variado, espaço ativo maior, geometria real,
   saltos pelo espaço, B3LYP e PBE0.
4. **O C₅₁H₂₅ reage consigo mesmo.** O xTB forma uma ligação σ C7–C26 (1,558 Å) no fiorde, o que
   apaga o magnetismo. A diferença grande que a DFT mostrava vinha de ter calculado essa outra
   molécula. A forma helicoidal sem a ligação existe, mas está 0,84 eV acima.
5. **Bioquímicas alternativas** acrescentam ~36% de civilizações às de carbono/água; o silício, só ~5%.

**O que ainda falta:** saber se o C₄₉H₂₁ já foi descrito na literatura. Ele é único dentro do
modelo, mas isso não prova que seja inédito. A PubChem está bloqueada neste ambiente.

## Camada A — surgimento

| cenário | chega à ciência | mediana após formação | idade do universo (mediana) |
|---|---|---|---|
| Rápido | 97,8% | 2,06 Gyr | 7,6 Gyr |
| Terra | 85,5% | 3,56 Gyr | 9,4 Gyr |
| Passos difíceis | 3,2% | 11,86 Gyr | 17,6 Gyr |

### Bioquímicas (cenário Terra, catálise adaptada)

| bioquímica | civilizações por civ. carbono/água | mediana até a ciência | T típica |
|---|---|---|---|
| Carbono/água | 1 | 3,7 Gyr | 308 K |
| Carbono/amônia | 0,141 | 5,8 Gyr | 213 K |
| Carbono/ácido sulfúrico | 0,120 | 3,5 Gyr | 397 K |
| Carbono/metano (Titã) | 0,047 | 11,1 Gyr | 136 K |
| Silício/água | 0,023 | 10,6 Gyr | 322 K |
| Silício/ácido sulfúrico | 0,023 | 9,4 Gyr | 430 K |
| Carbono/fluoreto de hidrogênio | 0,010 | 5,6 Gyr | 231 K |
| Silício/metano | 0,0014 | 12,6 Gyr | 142 K |
| Silício/nitrogênio líquido | 0,000007 | 8,6 Gyr | 71 K |

Com química fixa (Ea = 50 kJ/mol), as bioquímicas frias (amônia, metano, N₂) nunca chegam à ciência.
Energias de ligação, faixas líquidas e constantes dielétricas são dados reais; compatibilidade,
abundância dos solventes e o efeito da falta de fogo são hipóteses tiradas da literatura.
Detalhes: `resultados/bioquimica/resumo.md`.

## Camada B — civilizações pesquisando

- **Culturas humanas:** a isolada chega ao melhor resultado 35% a 60% mais tarde que as que trocam
  conhecimento, em todos os laboratórios. Os Empiristas fazem mais descobertas inéditas.
- **Bioquímicas (materiais):** carbono/ácido sulfúrico (pensamento 1,37×) chega ao recorde no ano 65;
  carbono/metano (0,48×, sem fogo) só no ano 116.
- **Química com filtro de planaridade:** a cultura isolada ficou bem atrás (score 4,98 contra 6,12).

## Química: verificação dos candidatos

Pipeline `scripts/pipeline_quimica.py`, sempre a partir do grafo exato exportado pelo simulador.

### Validação do método (Hubbard-CAS, t = 2,7 eV, U = 1,2·t)

| molécula | conhecido | calculado |
|---|---|---|
| Fenalenil C₁₃H₉ | S = 1/2 | S = 1/2 |
| [3]Triangulena C₂₂H₁₂ | S = 1 (Pavliček et al. 2017) | S = 1 |
| Cálice de Clar C₃₈H₁₈ | S = 0; tripleto +23 meV (Mishra et al. 2020) | S = 0; tripleto +23,5 meV |

O limiar de planaridade (H–H ≥ 1,5 Å) também foi validado: fenantreno e cálice de Clar passam
(1,75 Å); [4]- e [5]heliceno são rejeitados (0,57 e 0,76 Å).

### Candidatos

| | C₅₁H₂₅ (1ª busca) | C₄₉H₂₁ (busca plana) |
|---|---|---|
| origem | Empiristas, universo 1, ano 63 | Evolucionistas, universo 1, ano 218; final em 86/96 universos (3 sementes) |
| InChIKey | `LSGYJUYIXKPRHI-UHFFFAOYSA-N` | `SZADRCVFNYKXMY-UHFFFAOYSA-N` |
| Hückel | η = 3 (1 por desbalanço, 2 ocultos) | η = 3 (1 por desbalanço, 2 ocultos) |
| Hubbard | S = 1/2; S = 3/2 a +13,9 meV | S = 1/2; S = 3/2 a +21,3 meV |
| momento local máximo | 0,35 μ_B | 0,32 μ_B |
| xTB: planaridade | torcido, 2,1 Å | plano, 0,03 Å |
| xTB: spin 3/2 − spin 1/2 | +0,55 eV | +0,10 eV |
| xTB: frequências | nenhuma imaginária | nenhuma imaginária (menor: 24,7 cm⁻¹) |
| DFT: estado fundamental | spin 1/2 (simetria quebrada) | spin 1/2 (simetria quebrada) |
| DFT: spin 3/2 − spin 1/2 | +334 meV (projetado +440) | +50 meV (projetado +82) |
| DFT calibrada pelo cálice de Clar (÷3,8) | não se aplica: calculada na forma ciclizada | **~21,5 meV: bate** com o Hubbard |

Registros completos, com anéis, coordenadas, grafo, genoma, semente e comando, estão em
`data/candidatos/`. A busca química reproduz bit a bit com a mesma semente.

### DFT: correção, validação e resultado

A primeira tentativa convergiu para soluções erradas (camada quase fechada e configurações excitadas).
Correção: **chute inicial com a densidade de spin do Hubbard**, mais análise de estabilidade quando há
sinal de problema.

**Validação no cálice de Clar:** o estado fundamental (singleto de simetria quebrada) saiu correto,
com soluções estáveis, mas a distância até o tripleto foi 87,5 meV (projeção de Yamaguchi) contra
23 meV medidos. O UB3LYP/6-31G superestima ~3,8×, um desvio comum de funcionais híbridos nessas
moléculas. Por isso o número quantitativo de referência continua sendo o Hubbard-CAS, e a DFT confirma
estado fundamental, geometria e distribuição de spin.

**C₄₉H₂₁:** a razão DFT/Hubbard é 3,8×, a mesma do cálice de Clar. Três métodos diferentes concordam:
spin 1/2 fundamental, spin 3/2 a ~21 meV (~250 K). A densidade de spin mostra dois domínios
antiparalelos, os dois triângulos da molécula, como no cálice de Clar
(`resultados/candidatos/C49H21_eta3_plano/densidade_spin_dft.png`).

**C₅₁H₂₅:** ver a seção de robustez abaixo: a molécula cicliza, e a DFT foi feita na forma ciclizada.

Limites: base pequena (6-31G) e estabilidade verificada só no cálice de Clar. Um
cálculo publicável pediria base maior (def2-TZVP), outros funcionais e métodos multirreferência
(CASSCF/NEVPT2).

### Robustez (4ª rodada)

**Gap de spin (meV) em cada variação do Hubbard-CAS:**

| variação | cálice de Clar | C₄₉H₂₁ | razão |
|---|---|---|---|
| U/t = 1,0 | 16,7 | 15,0 | 0,90 |
| U/t = 1,2 (padrão) | 23,5 | 21,3 | 0,91 |
| U/t = 1,4 | 31,3 | 28,5 | 0,91 |
| U/t = 1,6 | 39,9 | 36,5 | 0,92 |
| espaço ativo maior (η+6 orbitais) | 19,5 | 18,9 | 0,97 |
| geometria xTB, saltos Slater–Koster só entre vizinhos (t·cos φ) | 18,7 | 14,8 | 0,79 |
| geometria xTB + saltos pelo espaço até 4 Å | 20,9 | 18,3 | 0,87 |

**DFT com dois funcionais (gap projetado de Yamaguchi):**

| | cálice de Clar | C₄₉H₂₁ | razão |
|---|---|---|---|
| UB3LYP/6-31G | 87,5 | 81,7 | 0,93 |
| UPBE0/6-31G | 140,7 | 126,9 | 0,90 |

O valor absoluto da DFT depende do funcional (mais troca exata, gap maior); a razão não depende.
A concordância exata de 23,5 meV no cálice de Clar com o CAS pequeno era em parte sorte: com o espaço
ativo maior ela vai a 19,5 meV. A razão entre as moléculas, que é o que importa, fica estável.

**C₅₁H₂₅, o que realmente aconteceu:** o xTB formou uma ligação σ C7–C26 de 1,558 Å, com os dois
carbonos virando sp³. É uma recombinação dos dois centros de maior densidade de spin, forçados a ficar
próximos no fiorde. A forma helicoidal sem a ligação (C7–C26 a 2,82 Å) é um mínimo, mas está 0,84 eV
acima da ciclizada. Nela, a torção sozinha (t·cos φ) reduz o gap de 13,9 para 9,4 meV, e o salto
direto entre as pontas da hélice o eleva para ~50 meV. A sugestão de explicar a discrepância só com
t·cos φ não se confirma.

**Isômeros (3 sementes, 96 universos):**

| molécula | InChIKey | resultado final (vezes) | universos | gap Hubbard |
|---|---|---|---|---|
| C₄₉H₂₁ | SZADRCVFNYKXMY-UHFFFAOYSA-N | 265 | 86 | 21,3 meV |
| C₃₈H₁₈ (cálice de Clar) | IQSWLUPJXFHUKA-UHFFFAOYSA-N | 6 | 4 | 23,5 meV |
| C₅₅H₂₅ | ONDYYYOZKMTRRY-UHFFFAOYSA-N | 4 | 2 | 7,5 meV |
| C₅₀H₂₀ | JOIPVWSJTCAJRK-UHFFFAOYSA-N | 2 | 2 | 24,1 meV |
| C₅₃H₂₃ | CRCBVLJGTUJWTH-UHFFFAOYSA-N | 1 | 1 | 16,6 meV |
| C₅₃H₂₅ | HUEHVUCUZKSIOZ-UHFFFAOYSA-N | 1 | 1 | 6,4 meV |

**Enumeração exaustiva** (`simu enumerar --aneis 14`, 13 min; `resultados/enumeracao/`):

| anéis | benzenoides (livres) | planos | com modos ocultos | com N_A = N_B | classe do C₄₉H₂₁ |
|---|---|---|---|---|---|
| ≤ 10 | 38 959 | 6 287 | 0 | 0 | 0 |
| 11 | 143 552 | 14 512 | 1 (cálice de Clar) | 1 | 0 |
| 12 | 683 101 | 48 866 | 4 | 4 | 0 |
| 13 | 3 274 826 | 165 789 | 37 | 37 | 0 |
| 14 | 15 796 897 | 565 320 | 237 | 236 | **1 (C₄₉H₂₁)** |

Validação: as contagens de formas fixas e livres batem exatamente com a OEIS (A001207 e A000228) em
todos os tamanhos, e o único caso de 11 anéis é o cálice de Clar, o menor da família segundo a
literatura. Um viés declarado: a pontuação da busca favorece moléculas compactas, e por isso a
enumeração era o teste decisivo.

### Literatura

A classe é conhecida: frustração topológica, a família do cálice de Clar. Um cálice de Clar estendido
(C₇₆H₂₆) foi sintetizado em 2025, e esses benzenoides são menos de 0,1% dos que têm até 14 anéis.
**A topologia exata não pôde ser verificada**: a busca web pelas InChIKeys não retornou nada, mas
PubChem e RSC estão bloqueados pela rede do ambiente. A checagem definitiva precisa das InChIKeys
acima na PubChem ou no CAS SciFinder.

## Materiais e nuclear

- **Materiais:** as bandas planas isoladas têm gap provado matematicamente (0 violações em 600 mil
  pontos k). A rodada completa caiu de 531 s para 39 s. Os mecanismos encontrados já são conhecidos.
- **Nuclear:** qualitativo. A ilha de estabilidade aparece em Mc-300 e Fl-290, mas as meias-vidas
  erram por várias ordens de grandeza.

## Próximos passos

1. Liberar `pubchem.ncbi.nlm.nih.gov` na rede do ambiente e checar a InChIKey do C₄₉H₂₁.
2. Estender a enumeração a 15–16 anéis e incluir não planos, para mapear toda a família.
3. DFT com base maior (def2-TZVP) e um método multirreferência (CASSCF/NEVPT2) no C₄₉H₂₁.
4. Estabilidade química: índices de reatividade dos sítios de spin (como no C₅₁H₂₅) e grupos
   protetores (mesitila), como se faz na síntese real de nanografenos radicalares.

## Onde estão os dados

| caminho | conteúdo |
|---|---|
| `resultados/*/resumo.md`, `*.csv`, `catalogo.json` | cada simulação completa |
| `resultados/bioquimica/` | surgimento por bioquímica |
| `resultados/*_bioquimicas/` | civilizações de bioquímicas diferentes pesquisando |
| `resultados/candidatos/` | Hubbard, geometrias xTB, DFT (B3LYP e PBE0), densidades de spin, robustez e isômeros |
| `resultados/enumeracao/` | enumeração exaustiva e lista de todos os planos com modos ocultos |
| `data/candidatos/` | estrutura exata dos candidatos |
| `resultados/relatorio.html` | relatório visual da 2ª rodada |
