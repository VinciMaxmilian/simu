# Relatório do simu — 3ª rodada

24 set 2026 · branch `claude/determined-tesla-nkx3fz` · também na `main`

## Resumo

1. **O magnetismo dos candidatos sobrevive à interação elétron-elétron.** No modelo de Hubbard
   (validado contra o experimento do cálice de Clar: 23,5 meV calculado × 23 meV medido), os dois
   candidatos têm estado fundamental de spin 1/2 e um estado de spin 3/2 logo acima.
2. **O C₅₁H₂₅ original não é plano.** Ele contém uma enseada de [4]heliceno e um fiorde de
   [5]heliceno (H–H de 0,57 e 0,76 Å na geometria plana). O xTB torce a molécula 2,1 Å. O laboratório
   de química agora rejeita estruturas não planas.
3. **Com o filtro de planaridade surgiu um candidato melhor, o C₄₉H₂₁**: plano (0,03 Å no xTB),
   mínimo estável (sem frequências imaginárias) e com acoplamento de spin mais forte (21,3 meV),
   encontrado em 28 dos 32 universos.
4. **A DFT ainda não é confiável.** Os cálculos convergiram para soluções erradas (detalhes abaixo). É
   o próximo passo do trabalho.
5. **Bioquímicas alternativas** acrescentam ~36% de civilizações às de carbono/água; o silício, só ~5%.

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
| origem | Empiristas, universo 1, ano 63 | Evolucionistas, universo 1, ano 218 |
| InChIKey | `LSGYJUYIXKPRHI-UHFFFAOYSA-N` | `SZADRCVFNYKXMY-UHFFFAOYSA-N` |
| Hückel | η = 3 (1 por desbalanço, 2 ocultos) | η = 3 (1 por desbalanço, 2 ocultos) |
| Hubbard | S = 1/2; S = 3/2 a +13,9 meV | S = 1/2; S = 3/2 a +21,3 meV |
| momento local máximo | 0,35 μ_B | 0,32 μ_B |
| xTB: planaridade | torcido, 2,1 Å | plano, 0,03 Å |
| xTB: spin 3/2 − spin 1/2 | +0,55 eV | +0,10 eV |
| xTB: frequências | nenhuma imaginária | nenhuma imaginária (menor: 24,7 cm⁻¹) |
| DFT (UB3LYP/6-31G) | **não confiável** | **não confiável** |

Registros completos, com anéis, coordenadas, grafo, genoma, semente e comando, estão em
`data/candidatos/`. A busca química reproduz bit a bit com a mesma semente.

### Por que a DFT ainda não vale

- **C₄₉H₂₁:** o spin 1/2 convergiu para uma solução quase de camada fechada (⟨S²⟩ = 0,81), e não
  para o estado de simetria quebrada que o Hubbard indica como fundamental. Com isso, o quarteto sai
  657 meV "abaixo", o que é artefato.
- **Quartetos dos dois candidatos:** o PySCF avisou "HOMO ≥ LUMO", ou seja, convergiram para uma
  configuração excitada. Isso infla o gap do C₅₁H₂₅ (461 meV contra 13,9 meV no Hubbard).

Correção em andamento:
1. análise de estabilidade da função de onda, reiniciando o cálculo a partir das instabilidades;
2. chute de simetria quebrada guiado pela densidade de spin do Hubbard;
3. validar a DFT no cálice de Clar (J experimental de 23 meV) antes de confiar nos candidatos.

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

1. Corrigir e validar a DFT (acima).
2. Levar os candidatos planos com maior acoplamento de spin a um critério de busca: avaliar os
   melhores de cada universo com Hubbard-CAS dentro da simulação.
3. Liberar `pubchem.ncbi.nlm.nih.gov` na rede do ambiente para a checagem de literatura.

## Onde estão os dados

| caminho | conteúdo |
|---|---|
| `resultados/*/resumo.md`, `*.csv`, `catalogo.json` | cada simulação completa |
| `resultados/bioquimica/` | surgimento por bioquímica |
| `resultados/*_bioquimicas/` | civilizações de bioquímicas diferentes pesquisando |
| `resultados/candidatos/` | Hubbard, geometrias xTB, DFT, densidades de spin e logs |
| `data/candidatos/` | estrutura exata dos candidatos |
| `resultados/relatorio.html` | relatório visual da 2ª rodada |
