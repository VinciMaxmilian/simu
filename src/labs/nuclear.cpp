// Fase 3 — Física nuclear: quais núcleos (Z prótons, N nêutrons) existem e quanto vivem?
//
// Física (semi-empírica, derivada da QCD + eletromagnetismo do Modelo Padrão só de forma efetiva):
//   • energia de ligação: gota líquida (Bethe–Weizsäcker) + correção de camadas de Myers–Swiatecki
//     com números mágicos (os de Z=114/164 e N=184/258 são PREVISÕES teóricas, entradas do modelo);
//   • decaimento β: regra de Sargent (T ∝ Q⁻⁵);
//   • decaimento α: fórmula de Viola–Seaborg (parâmetros de Sobiczewski 1989);
//   • fissão espontânea: barreira da gota líquida + correção de camadas, log10 T = 2.30·B_f − 14.9.
// Calibração (scripts/calibrar_nuclear.py): erro rms 78 keV em B/A, 0.83 MeV em Qα, 4 décadas em T_fissão.
// Como Qα entra exponencialmente na meia-vida α, as meias-vidas têm erro de várias ordens de grandeza:
// este laboratório é QUALITATIVO (onde fica a ilha, não quanto tempo exato ela vive).
// Objetivo: o elemento mais pesado com meia-vida ≥ 60 s ("ilha de estabilidade").
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "labs/labs.hpp"

namespace {

constexpr int kMaxZ = 140;
constexpr int kMaxN = 230;
constexpr double kMinLog10T = -21.0;  // "não existe"

const char* kSymbols[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne", "Na", "Mg", "Al", "Si", "P",  "S",
    "Cl", "Ar", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho",
    "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
    "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md",
    "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};

std::string symbol(int z) { return z <= 118 ? kSymbols[z] : "E" + std::to_string(z); }

// Função de camadas de Myers–Swiatecki: 0 nos números mágicos, máxima no meio da camada.
double shell_f(double x, const std::vector<double>& magic) {
  double prev = 0;
  for (double m : magic) {
    if (x <= m) {
      const double q = 0.6 * (std::pow(m, 5.0 / 3) - std::pow(prev, 5.0 / 3)) / (m - prev);
      return q * (x - prev) - 0.6 * (std::pow(x, 5.0 / 3) - std::pow(prev, 5.0 / 3));
    }
    prev = m;
  }
  return 0;
}

const std::vector<double> kMagicZ = {2, 8, 20, 28, 50, 82, 114, 164};
const std::vector<double> kMagicN = {2, 8, 20, 28, 50, 82, 126, 184, 258};

double shell_correction(int z, int n) {
  const double a = z + n;
  const double s = 6.11388 * (shell_f(z, kMagicZ) + shell_f(n, kMagicN)) / std::pow(a / 2, 2.0 / 3) -
                   1.81224 * std::cbrt(a);
  return std::min(s, 0.0);  // calibração escolheu k_deformação = 0: S > 0 é anulado por deformação
}

double binding_ld(int z, int n) {
  const double a = z + n;
  // Coeficientes ajustados por scripts/calibrar_nuclear.py (dados em data/nuclear/calibracao.csv).
  const double i2 = static_cast<double>(n - z) * (n - z);
  double b = 15.67692 * a - 18.09073 * std::pow(a, 2.0 / 3) - 0.70311 * z * (z - 1) / std::cbrt(a) -
             27.74643 * i2 / a + 28.35154 * i2 / std::pow(a, 4.0 / 3);
  if (z % 2 == 0 && n % 2 == 0) b += 12.0 / std::sqrt(a);
  if (z % 2 == 1 && n % 2 == 1) b -= 12.0 / std::sqrt(a);
  return b;
}

// Energia de ligação total (MeV). Núcleos sem sentido físico recebem −∞.
double binding(int z, int n) {
  if (z < 1 || n < 0 || z + n < 2) return -1e9;
  return binding_ld(z, n) - shell_correction(z, n);
}

struct Decay {
  double log10_t = kMinLog10T;  // meia-vida total (s)
  std::string mode;
  double t_alpha = 99, t_beta = 99, t_sf = 99;  // log10 de cada canal
};

Decay half_life(int z, int n) {
  Decay d;
  const double a = z + n;
  const double b = binding(z, n);
  if (b <= 0 || b - binding(z, n - 1) <= 0 || b - binding(z - 1, n) <= 0) {
    d.mode = "instável a emissão de núcleons";
    return d;
  }
  auto add = [](double acc, double log10_t) { return acc + std::pow(10.0, -log10_t); };
  double rate = 0;

  // β⁻ e β⁺/captura eletrônica (Q com massas atômicas: m_n − m_H = 0.782 MeV).
  const double q_minus = binding(z + 1, n - 1) - b + 0.782;
  const double q_plus = binding(z - 1, n + 1) - b - 0.782;
  const double q_beta = std::max(q_minus, q_plus);
  if (q_beta > 0) {
    d.t_beta = 4.0 - 5.0 * std::log10(q_beta);
    rate = add(rate, d.t_beta);
  }
  // α (Viola–Seaborg), com impedimento para núcleons ímpares.
  if (z > 2 && n > 2) {
    const double q_alpha = binding(z - 2, n - 2) + 28.296 - b;
    if (q_alpha > 0.5) {
      d.t_alpha = (1.66175 * z - 8.5166) / std::sqrt(q_alpha) - 0.20228 * z - 33.9069 + (z % 2) * 0.772 +
                  (n % 2) * 1.066;
      if (d.t_alpha < 40) rate = add(rate, d.t_alpha);
    }
  }
  // Fissão espontânea: barreira da gota líquida (Cohen–Swiatecki aprox.) + camadas do estado fundamental.
  if (z > 80) {
    const double es = 18.09073 * std::pow(a, 2.0 / 3), ec = 0.70311 * z * z / std::cbrt(a);
    const double x = ec / (2 * es);
    double bf_ld = 0;
    if (x < 2.0 / 3) bf_ld = 0.38 * (0.75 - x) * es;
    else if (x < 1) bf_ld = 0.83 * std::pow(1 - x, 3) * es;
    const double bf = bf_ld - shell_correction(z, n);
    d.t_sf = std::max(kMinLog10T, 2.3009 * bf - 14.945);
    if (d.t_sf < 40) rate = add(rate, d.t_sf);
  }
  if (rate == 0) {
    d.log10_t = 99;
    d.mode = "estável";
    return d;
  }
  d.log10_t = std::max(kMinLog10T, -std::log10(rate));
  const double m = std::min({d.t_alpha, d.t_beta, d.t_sf});
  d.mode = m == d.t_alpha ? "α" : m == d.t_beta ? "β" : "fissão";
  return d;
}

std::string format_time(double log10_t) {
  if (log10_t >= 99) return "estável";
  const double t = std::pow(10.0, log10_t);
  char buf[64];
  if (t < 1) std::snprintf(buf, sizeof buf, "%.1e s", t);
  else if (t < 3600) std::snprintf(buf, sizeof buf, "%.3g s", t);
  else if (t < 3.156e7) std::snprintf(buf, sizeof buf, "%.3g h", t / 3600);
  else std::snprintf(buf, sizeof buf, "%.3g anos", t / 3.156e7);
  return buf;
}

class NuclearLab : public Lab {
 public:
  std::string id() const override { return "nuclear"; }
  std::string title() const override { return "Física nuclear: a ilha de estabilidade (modelo semi-empírico)"; }
  std::string objective() const override {
    return "achar o elemento mais pesado que vive ≥ 60 s; score = Z − 5 por década abaixo de 60 s";
  }

  Genome random(Rng& rng) const override {
    const int z = 1 + rng.below(kMaxZ);
    const int n = std::clamp(static_cast<int>(z * rng.uniform(0.9, 1.7)), 0, kMaxN);
    return {static_cast<uint8_t>(z), static_cast<uint8_t>(n)};
  }
  Genome mutate(const Genome& g, Rng& rng) const override {
    Genome m = g;
    const int which = rng.below(3);
    const int step = (1 + rng.below(3)) * (rng.chance(0.5) ? 1 : -1);
    if (which != 1) m[0] = static_cast<uint8_t>(std::clamp(m[0] + step, 1, kMaxZ));
    if (which != 0) m[1] = static_cast<uint8_t>(std::clamp(m[1] + step, 0, kMaxN));
    return m;
  }
  Genome crossover(const Genome& a, const Genome& b, Rng& rng) const override {
    return rng.chance(0.5) ? Genome{a[0], b[1]} : Genome{b[0], a[1]};
  }

  Evaluation evaluate(const Genome& g) const override {
    Evaluation e;
    const int z = g[0], n = g[1];
    if (z < 1 || z > kMaxZ || n > kMaxN || z + n < 2) return e;
    const Decay d = half_life(z, n);
    if (d.log10_t <= kMinLog10T) return e;
    e.valid = true;
    e.score = z - 5.0 * std::max(0.0, std::log10(60.0) - d.log10_t);
    e.label = d.log10_t >= std::log10(60.0) ? "Z=" + std::to_string(z) + " " + symbol(z) + " (vida longa)"
                                            : "nucleos efemeros";
    return e;
  }

  std::string describe(const Genome& g) const override {
    const int z = g[0], n = g[1];
    const Decay d = half_life(z, n);
    char buf[256];
    std::snprintf(buf, sizeof buf, "%s-%d (Z=%d, N=%d): B/A=%.3f MeV, meia-vida %s (%s)", symbol(z).c_str(), z + n, z,
                  n, binding(z, n) / (z + n), format_time(d.log10_t).c_str(), d.mode.c_str());
    return buf;
  }

  std::vector<std::string> validate() const override {
    std::vector<std::string> out;
    out.push_back("| núcleo | B/A real (MeV) | B/A modelo | meia-vida real | meia-vida modelo |");
    out.push_back("|---|---|---|---|---|");
    struct Ref { int z, n; double ba; const char* t; };
    const Ref refs[] = {{2, 2, 7.074, "estável"},     {6, 6, 7.680, "estável"},     {8, 8, 7.976, "estável"},
                        {20, 20, 8.551, "estável"},   {26, 30, 8.790, "estável"},   {28, 34, 8.795, "estável"},
                        {50, 70, 8.505, "estável"},   {82, 126, 7.867, "estável"},  {84, 128, 7.810, "0.3 μs"},
                        {84, 126, 7.834, "138 dias"}, {90, 142, 7.615, "1.4e10 anos"},
                        {92, 146, 7.570, "4.5e9 anos"}, {94, 145, 7.560, "2.4e4 anos"}, {100, 156, 7.4, "2.6 h"},
                        {114, 175, 7.0, "~2 s"},      {118, 176, 6.9, "0.7 ms"}};
    for (const auto& r : refs) {
      const Decay d = half_life(r.z, r.n);
      char buf[256];
      std::snprintf(buf, sizeof buf, "| %s-%d | %.3f | %.3f | %s | %s (%s) |", symbol(r.z).c_str(), r.z + r.n, r.ba,
                    binding(r.z, r.n) / (r.z + r.n), r.t, format_time(d.log10_t).c_str(), d.mode.c_str());
      out.push_back(buf);
    }
    // Pico de B/A e elemento estável mais pesado segundo o modelo.
    double best_ba = 0;
    int bz = 0, bn = 0, heaviest_stable = 0;
    for (int z = 1; z <= kMaxZ; ++z)
      for (int n = 0; n <= kMaxN; ++n) {
        if (z + n < 2) continue;
        const double ba = binding(z, n) / (z + n);
        if (ba > best_ba) best_ba = ba, bz = z, bn = n;
        if (z > heaviest_stable && half_life(z, n).log10_t >= 99) heaviest_stable = z;
      }
    out.push_back("");
    out.push_back("Pico de B/A no modelo: " + symbol(bz) + "-" + std::to_string(bz + bn) +
                  " (real: Ni-62 8.795 MeV, Fe-56 8.790 MeV).");
    out.push_back("Elemento estável mais pesado no modelo: Z=" + std::to_string(heaviest_stable) + " " +
                  symbol(heaviest_stable) + " (real: Pb, Z=82; Bi-209 decai em 2e19 anos).");
    return out;
  }
};

}  // namespace

std::unique_ptr<Lab> make_nuclear_lab() { return std::make_unique<NuclearLab>(); }
