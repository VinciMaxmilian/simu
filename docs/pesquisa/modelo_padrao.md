# Modelo Padrão: compêndio de dados e parâmetros

Fontes: PDG 2024/2025 (Review of Particle Physics), CODATA 2022, NuFIT 6.0 (2024), Fermilab g-2 (2025).
⚠️ O acesso ao pdg.lbl.gov está bloqueado neste ambiente. Os valores marcados com (≈) vieram de memória e
de buscas e precisam ser conferidos na implementação.

## A1. Conteúdo de partículas (17 partículas fundamentais)

**Férmions (spin 1/2), 3 gerações**

| Partícula | Símbolo | Carga (e) | Cor | Massa | Descoberta |
|---|---|---|---|---|---|
| up | u | +2/3 | sim | 2.16 MeV (+0.49 −0.26) (MS̄, 2 GeV) | 1968 (SLAC, espalhamento inelástico profundo) |
| down | d | −1/3 | sim | 4.70 MeV (≈ ±0.07) | 1968 (SLAC) |
| charm | c | +2/3 | sim | 1.2730 ± 0.0046 GeV (MS̄) | 1974 (J/ψ, SLAC+BNL) |
| strange | s | −1/3 | sim | 93.5 ± 0.8 MeV | 1947 (kaons) / 1968 |
| top | t | +2/3 | sim | 172.56 ± 0.31 GeV (direta, PDG 2025) | 1995 (Tevatron CDF/D0) |
| bottom | b | −1/3 | sim | 4.183 ± 0.007 GeV (MS̄) | 1977 (Υ, Fermilab) |
| elétron | e | −1 | não | 0.51099895000(15) MeV | 1897 (Thomson) |
| múon | μ | −1 | não | 105.6583755(23) MeV | 1936 (Anderson & Neddermeyer) |
| tau | τ | −1 | não | 1776.93 ± 0.09 MeV | 1975 (Perl, SLAC) |
| neutrino-e | νe | 0 | não | < 0.45 eV (KATRIN, 90% CL) | 1956 (Cowan–Reines) |
| neutrino-μ | νμ | 0 | não | Σmν < 0.12 eV (Planck); < ~0.07 (DESI, ≈) | 1962 (BNL) |
| neutrino-τ | ντ | 0 | não | — | 2000 (DONUT) |

**Bósons**

| Partícula | Spin | Carga | Massa | Largura | Descoberta |
|---|---|---|---|---|---|
| fóton γ | 1 | 0 | 0 (< 1e-18 eV) | — | 1905 (Einstein) / 1923 (Compton) |
| glúon g (8) | 1 | 0 | 0 | — | 1979 (PETRA, eventos de 3 jatos) |
| W± | 1 | ±1 | 80.3692 ± 0.0133 GeV (PDG 2024); ~80.3625 ± 0.0077 sem CDF (≈ 2025) | 2.085 ± 0.042 GeV | 1983 (UA1/UA2, CERN) |
| Z⁰ | 1 | 0 | 91.1880 ± 0.0020 GeV | 2.4955 ± 0.0023 GeV | 1983 (CERN) |
| Higgs H | 0 | 0 | 125.20 ± 0.11 GeV | ~3.7 MeV (+1.9 −1.4) medida; 4.1 MeV no MP | 2012 (ATLAS/CMS) |

**Números quânticos** (convenção Q = T₃ + Y):
Q_L (1/6, dubleto), u_R (2/3), d_R (−1/3), L_L (−1/2), e_R (−1), H (1/2).
Número bariônico: 1/3 por quark. Número leptônico: 1 por lépton. Toda partícula tem antipartícula (o fóton, o Z, o H e os glúons neutros são as próprias).

**Grupo de gauge:** SU(3)_c × SU(2)_L × U(1)_Y → SU(3)_c × U(1)_em após a quebra de simetria.

## A2. Os 19 parâmetros livres (+7 a 9 com neutrinos massivos)

| # | Parâmetro | Valor |
|---|---|---|
| 1–6 | massas dos quarks | tabela acima |
| 7–9 | massas dos léptons carregados | tabela acima |
| 10–12 | ângulos CKM (Wolfenstein) | λ = 0.22501 ± 0.00068; A = 0.826 (+0.016 −0.015); ρ̄ = 0.1591 ± 0.0094 |
| 13 | fase CKM | η̄ = 0.3523 (+0.0073 −0.0071); J (Jarlskog) ≈ 3.08e-5 |
| 14 | α_s(M_Z) | 0.1180 ± 0.0009 |
| 15 | α (em, Q=0) | 1/137.035999177(21) (CODATA 2022) |
| 16 | sin²θ_W (MS̄, M_Z) | 0.23129 ± 0.00004 (≈) ; on-shell ≈ 0.2232 |
| 17 | θ_QCD | < ~1e-10 (EDM do nêutron \|d_n\| < 1.8e-26 e·cm) |
| 18 | massa do Higgs | 125.20 ± 0.11 GeV |
| 19 | vev do Higgs v | 246.22 GeV (= (√2 G_F)^(-1/2)) |

Derivados úteis: G_F = 1.1663788(6)e-5 GeV⁻²; α(M_Z) ≈ 1/127.95; g ≈ 0.652, g' ≈ 0.357, g_s(M_Z) ≈ 1.22;
λ_Higgs ≈ 0.13; Yukawa do top ≈ 0.99.

**Matriz CKM (magnitudes, PDG ≈):**
```
|Vud| 0.97367(32)  |Vus| 0.22431(85)  |Vub| 3.82(20)e-3
|Vcd| 0.221(4)     |Vcs| 0.975(6)     |Vcb| 41.0(1.4)e-3
|Vtd| 8.6(2)e-3    |Vts| 41.5(9)e-3   |Vtb| 1.010(27)
```

**Neutrinos (NuFIT 6.0, ordenamento normal):** θ12 = 33.68° (sin² = 0.307); θ23 = 43.3° (sin² ≈ 0.470,
octante ambíguo); θ13 = 8.56° (sin² ≈ 0.0221); δ_CP = 212° (+26 −41); Δm²21 ≈ 7.49e-5 eV²;
Δm²31 ≈ +2.513e-3 eV² (≈). O ordenamento (normal ou invertido) segue desconhecido. As fases de Majorana e a natureza Dirac/Majorana também.

## A3. Hádrons e vidas médias de referência

próton 938.27208816 MeV, estável (τ > 2.4e34 anos, Super-K), raio 0.8409 fm · nêutron 939.56542052 MeV,
τ = 878.4 ± 0.5 s · π± 139.57 MeV, τ = 26.03 ns · π⁰ 134.98 MeV, τ = 8.4e-17 s · K± 493.68 MeV ·
múon τ = 2.1969811(22) μs · tau τ = 290.3 fs · top τ ≈ 5e-25 s (não hadroniza).

## A4. Constantes físicas (SI 2019, exatas salvo indicação)

c = 299 792 458 m/s · h = 6.62607015e-34 J·s · e = 1.602176634e-19 C · k_B = 1.380649e-23 J/K ·
N_A = 6.02214076e23 · G = 6.67430(15)e-11 (medida) · ħc = 197.3269804 MeV·fm · M_Planck = 1.22089e19 GeV.

## A5. Problemas em aberto (nós "ainda não descobertos" no simulador)

matéria escura (~27% do universo) · energia escura (~68%) · origem da massa dos neutrinos · assimetria
matéria-antimatéria (η_B ≈ 6.1e-10) · gravidade quântica · problema da hierarquia · CP forte (áxion?) ·
por que 3 gerações · estabilidade do vácuo (m_H e m_t indicam metaestabilidade) · unificação.
Anomalias: o g-2 do múon foi resolvido em 2025 (exp. 116592070.5(14.8)e-11 × teoria WP25 116592033(62)e-11,
em acordo). A massa do W da CDF (80.4335) é uma medida isolada, incompatível com a média.

## A6. Linha do tempo histórica (calibração da simulação)

Os pares "previsão → confirmação" medem a defasagem entre teoria e experimento:

| Nó | Previsto | Confirmado | Pré-requisito tecnológico |
|---|---|---|---|
| elétron | — | 1897 | tubo de raios catódicos, vácuo |
| núcleo | — | 1911 (Rutherford) | fontes radioativas |
| fóton | 1905 | 1923 | espectroscopia |
| próton | — | 1919 | espalhamento α |
| pósitron | 1928 (Dirac) | 1932 | câmara de nuvens + raios cósmicos |
| neutrino | 1930 (Pauli) | 1956 | reator nuclear |
| nêutron | 1920 (Rutherford) | 1932 (Chadwick) | — |
| méson (píon) | 1935 (Yukawa) | 1947 | emulsões fotográficas |
| múon | — (inesperado) | 1936 | raios cósmicos |
| quarks | 1964 (Gell-Mann/Zweig) | 1968 | acelerador linear de 20 GeV |
| charm | 1970 (GIM) | 1974 | colisores e⁺e⁻ |
| 3ª geração | 1973 (Kobayashi–Maskawa) | 1975–1995 | — |
| glúon / QCD | 1973 (liberdade assintótica) | 1979 | colisor PETRA |
| W/Z | 1967–68 (Glashow–Weinberg–Salam) | 1983 | colisor p p̄ + resfriamento estocástico |
| top | 1973 | 1995 | Tevatron (~2 TeV) |
| oscilação de neutrinos | 1957–68 (Pontecorvo) | 1998 / 2001 | Super-K, SNO |
| Higgs | 1964 | 2012 | LHC (~8 TeV) |

Outros marcos: renormalizabilidade ('t Hooft, 1971), violação de CP (1964), violação de paridade (1956, Wu).

## Fontes

- Particle Data Group, *Review of Particle Physics* 2024/2025: https://pdg.lbl.gov (massa do top 172.56 ± 0.31 GeV, PDG 2025; Higgs 125.20 ± 0.11 GeV; W sem CDF ≈ 80.3625 GeV)
- CODATA 2022 (constantes fundamentais)
- NuFIT 6.0 (2024), arXiv:2410.05380 (oscilação de neutrinos)
- Fermilab Muon g-2, resultado final (junho de 2025) e Theory Initiative WP25

⚠️ O acesso direto ao pdg.lbl.gov estava bloqueado no ambiente em que isto foi escrito. Os valores marcados com (≈), e os campos `"verificar": true` em `data/modelo_padrao/*.json`, precisam ser conferidos na fonte.

