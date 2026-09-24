// Fase 1 — Materiais: cristais 2D em modelo tight-binding, à procura de BANDAS PLANAS.
//
// Por quê: numa banda plana a energia cinética dos elétrons some e as interações dominam —
// é a origem da supercondutividade e do magnetismo no grafeno bicamada torcido (2018) e em
// metais de rede kagome. Achar redes com bandas planas ISOLADAS (separadas das demais por um
// gap) é um problema de pesquisa ativo.
//
// Física: elétrons numa rede periódica (Bloch). Célula unitária com 1–4 sítios numa rede de
// Bravais quadrada; saltos entre vizinhos até (±1, ±1) com amplitude t ∈ {−1, +1};
// energias locais ε ∈ {−1, 0, +1}.
//
// Como a avaliação é rigorosa e rápida:
//   • triagem em 2 estágios: a largura de uma banda numa grade é sempre ≤ a largura real, então
//     largura > 0.05 numa grade 4×4 já PROVA que a banda não é plana (maioria das redes para aí);
//   • simetria de reversão temporal (saltos reais ⇒ H(−k) = H(k)*): metade da grade basta;
//   • gap CERTIFICADO: det(H(k) − E) é um polinômio trigonométrico de grau ≤ n em kx e ky
//     (saltos só até células vizinhas), determinado exatamente por (2n+1)² amostras. Com seus
//     coeficientes de Fourier, um branch-and-bound com limites de Taylor rigorosos prova que
//     det ≠ 0 em TODA a zona de Brillouin — i.e., nenhuma banda passa pela energia E. Provado isso
//     para dois níveis E1 < E2 dentro do gap aparente, o gap real é ≥ E2 − E1. Sem falsos isolados.
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <set>

#include "core/linalg.hpp"
#include "labs/labs.hpp"

namespace {

constexpr int kMaxSites = 4;
constexpr int kScreenGrid = 4;  // triagem: largura aqui é limite inferior da largura real
constexpr int kGrid = 9;        // candidatos: 9 = 2n+1 amostras por eixo; planura nesta grade implica
                                // det(H − E0) ≡ 0 (polinômio trigonométrico de grau ≤ 4 zerado em 9×9 nós)

struct Bond { int i, j, dx, dy; };

// Lista fixa de ligações candidatas: dentro da célula (i<j) e para 4 células vizinhas (todos os pares).
std::vector<Bond> candidate_bonds() {
  std::vector<Bond> out;
  for (int i = 0; i < kMaxSites; ++i)
    for (int j = i + 1; j < kMaxSites; ++j) out.push_back({i, j, 0, 0});
  const int disp[4][2] = {{1, 0}, {0, 1}, {1, 1}, {1, -1}};
  for (const auto& d : disp)
    for (int i = 0; i < kMaxSites; ++i)
      for (int j = 0; j < kMaxSites; ++j) out.push_back({i, j, d[0], d[1]});
  return out;
}

const std::vector<Bond>& bonds() {
  static const std::vector<Bond> b = candidate_bonds();
  return b;
}

// Layout do genoma: [n−1] [ε_0..ε_3] [ligação_0 .. ligação_B−1]
//   ε: 0 → 0, 1 → −1, 2 → +1   |   ligação: 0 → ausente, 1 → t=−1, 2 → t=+1
constexpr int kOnsite = 1;
constexpr int kBonds = 1 + kMaxSites;

struct Lattice {
  int n;
  double onsite[kMaxSites];
  std::vector<std::pair<Bond, double>> hops;
};

struct Bands {
  int n = 0;
  std::vector<std::vector<double>> e;  // [k][banda], só um representante de cada par (k, −k)
};

// f(kx, ky) = Σ c_ab e^{i(a·kx + b·ky)}, a, b ∈ [−deg, deg]; f real.
struct TrigPoly {
  int deg = 0;
  std::vector<std::complex<double>> c;  // [(a+deg)·(2deg+1) + (b+deg)]
  double l1 = 0, cubic = 0;  // Σ|c|  e  Σ|c|(a²+b²)^{3/2} (limite da 3ª derivada direcional)
  double sample_min = 0, sample_max = 0;  // valores nas amostras usadas na DFT

  // Valor, gradiente e hessiana num ponto (derivadas de e^{i(ax+by)}: fatores ia, ib).
  struct Local { double f, fx, fy, fxx, fxy, fyy; };
  Local eval(double x, double y) const {
    const int m = 2 * deg + 1;
    std::complex<double> ex[9], ey[9];
    const auto ux = std::polar(1.0, x), uy = std::polar(1.0, y);
    ex[deg] = ey[deg] = 1.0;
    for (int a = 1; a <= deg; ++a) {
      ex[deg + a] = ex[deg + a - 1] * ux, ex[deg - a] = std::conj(ex[deg + a]);
      ey[deg + a] = ey[deg + a - 1] * uy, ey[deg - a] = std::conj(ey[deg + a]);
    }
    // Somas por linha: R_a = Σ_b c_ab e^{iby}, S_a = Σ_b b·(…), Q_a = Σ_b b²·(…).
    Local out{0, 0, 0, 0, 0, 0};
    for (int a = -deg; a <= deg; ++a) {
      std::complex<double> r = 0, sb = 0, qb = 0;
      const std::complex<double>* ca = &c[(a + deg) * m];
      for (int b = -deg; b <= deg; ++b) {
        const auto term = ca[b + deg] * ey[b + deg];
        r += term, sb += term * static_cast<double>(b), qb += term * static_cast<double>(b * b);
      }
      const auto e = ex[a + deg];
      const auto re = r * e, se = sb * e, qe = qb * e;
      out.f += re.real();
      out.fx -= a * re.imag();
      out.fy -= se.imag();
      out.fxx -= a * a * re.real();
      out.fxy -= a * se.real();
      out.fyy -= qe.real();
    }
    return out;
  }
};

struct FlatInfo {
  int band = -1;
  double width = 1e9, gap = 0, grid_gap = 0, energy = 0, score = -1e9;
};

class MaterialsLab : public Lab {
 public:
  std::string id() const override { return "materiais"; }
  std::string title() const override { return "Materiais: redes 2D com bandas planas (tight-binding)"; }
  std::string objective() const override {
    return "achar a banda mais plana possível e, de preferência, isolada das outras por um gap de energia; "
           "score = −log10(largura + 1e−6) + 3·min(gap certificado, 1)";
  }

  Genome random(Rng& rng) const override {
    Genome g(kBonds + bonds().size(), 0);
    g[0] = static_cast<uint8_t>(rng.below(kMaxSites));
    for (int s = 0; s < kMaxSites; ++s) g[kOnsite + s] = rng.chance(0.15) ? 1 + rng.below(2) : 0;
    const double density = rng.uniform(0.03, 0.12);
    for (size_t b = 0; b < bonds().size(); ++b)
      if (rng.chance(density)) g[kBonds + b] = rng.chance(0.8) ? 1 : 2;
    return g;
  }

  Genome mutate(const Genome& g, Rng& rng) const override {
    Genome m = g;
    const int edits = 1 + rng.below(3);
    for (int e = 0; e < edits; ++e) {
      const double r = rng.uniform();
      if (r < 0.05) {
        m[0] = static_cast<uint8_t>(rng.below(kMaxSites));
      } else if (r < 0.20) {
        m[kOnsite + rng.below(kMaxSites)] = static_cast<uint8_t>(rng.below(3));
      } else {
        const size_t b = kBonds + rng.below(static_cast<int>(bonds().size()));
        m[b] = m[b] ? 0 : static_cast<uint8_t>(1 + rng.below(2));
      }
    }
    return m;
  }

  Genome crossover(const Genome& a, const Genome& b, Rng& rng) const override {
    Genome c(a.size());
    for (size_t i = 0; i < a.size(); ++i) c[i] = rng.chance(0.5) ? a[i] : b[i];
    return c;
  }

  Evaluation evaluate(const Genome& g) const override {
    Evaluation e;
    const Lattice lat = decode(g);
    if (!valid(lat)) return e;
    const FlatInfo f = analyze(lat, kScreenGrid);
    e.valid = true;
    e.score = f.score;
    e.label = classify(lat, f);
    return e;
  }

  // Só o que entra no hamiltoniano: sítios usados, suas energias e seus saltos.
  std::string canonical(const Genome& g) const override {
    const Lattice lat = decode(g);
    std::string key(1, static_cast<char>(lat.n));
    for (int i = 0; i < lat.n; ++i) key += static_cast<char>(lat.onsite[i] + 2);
    for (const auto& [b, t] : lat.hops)
      key += {static_cast<char>(b.i), static_cast<char>(b.j), static_cast<char>(b.dx + 2), static_cast<char>(b.dy + 2),
              static_cast<char>(t > 0)};
    return key;
  }

  std::string describe(const Genome& g) const override {
    const Lattice lat = decode(g);
    const FlatInfo f = analyze(lat, kScreenGrid, 24);
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "%d sítios/célula; banda %d de %d em E=%.3f|t|, largura=%.2e (grade 24×24), "
                  "gap certificado ≥ %.3f (aparente na grade: %.3f)",
                  lat.n, f.band + 1, lat.n, f.energy, f.width, f.gap, f.grid_gap);
    std::string s = buf;
    s += "\n    energias locais:";
    for (int i = 0; i < lat.n; ++i) s += " " + std::to_string(static_cast<int>(lat.onsite[i]));
    s += "\n    saltos (i→j, célula, t):";
    for (const auto& [b, t] : lat.hops) {
      std::snprintf(buf, sizeof buf, " %d→%d(%d,%d)%+d", b.i, b.j, b.dx, b.dy, static_cast<int>(t));
      s += buf;
    }
    return s;
  }

  std::vector<std::string> validate() const override {
    std::vector<std::string> out;
    auto show = [&](const char* name, const char* expected, const Lattice& lat, int grid = 24) {
      const FlatInfo f = analyze(lat, kScreenGrid, grid);
      const Bands b = compute_bands(lat, 24);
      char buf[320];
      std::snprintf(buf, sizeof buf,
                    "| %s | %s | banda mais plana: E=%.3f, largura=%.2e, gap certificado ≥ %.3f (aparente %.3f); "
                    "faixa total [%.3f, %.3f] → `%s` |",
                    name, expected, f.energy, f.width, f.gap, f.grid_gap, min_energy(b), max_energy(b),
                    classify(lat, f).c_str());
      out.push_back(buf);
    };
    out.push_back("| rede | esperado (teoria conhecida) | modelo |");
    out.push_back("|---|---|---|");
    Lattice square{1, {0}, {{{0, 0, 1, 0}, -1}, {{0, 0, 0, 1}, -1}}};
    show("quadrada", "uma banda, E ∈ [−4, 4], largura 8", square);
    Lattice lieb{3, {0, 0, 0}, {{{0, 1, 0, 0}, -1}, {{0, 2, 0, 0}, -1}, {{1, 0, 1, 0}, -1}, {{2, 0, 0, 1}, -1}}};
    show("Lieb (CuO₂ dos cupratos)", "banda plana em E=0 tocando as outras (gap 0) → não isolada", lieb);
    Lattice kagome{3, {0, 0, 0},
                   {{{0, 1, 0, 0}, -1}, {{0, 2, 0, 0}, -1}, {{1, 2, 0, 0}, -1},
                    {{1, 0, 1, 0}, -1}, {{2, 0, 0, 1}, -1}, {{1, 2, 1, -1}, -1}}};
    show("kagome", "banda plana em E=+2|t| tocando em Γ (gap 0) → não isolada", kagome);
    // Regressão: rede que a versão por amostragem classificava como "isolada" (toca fora da grade).
    Lattice fake{3, {0, 0, 0}, {{{0, 1, 1, 1}, -1}, {{1, 0, 1, 1}, 1}, {{1, 2, 1, -1}, 1}, {{2, 1, 1, -1}, 1}}};
    show("falso isolado da v1, analisado na grade 12×12 onde a v1 errava",
         "banda plana tocando (as bandas se encontram fora da grade 12×12)", fake, 12);
    Lattice isolated{3, {0, 0, 0}, {{{0, 1, 0, 1}, -1}, {{0, 2, 0, 1}, -1}, {{1, 0, 0, 1}, -1}, {{2, 0, 0, 1}, 1},
                                    {{0, 0, 1, 1}, -1}}};
    show("isolada achada pelas civilizações", "banda plana isolada (gap ≈ 1.236 na grade 24×24)", isolated);
    out.push_back("");
    out.push_back(stress_test());
    return out;
  }

 private:
  // Teste de estresse do certificado: redes aleatórias com gap certificado g; em milhares de pontos k
  // contínuos aleatórios, a banda escolhida não pode chegar a menos de g das vizinhas.
  std::string stress_test() const {
    Rng rng(2024);
    int certified = 0, violations = 0;
    long points = 0;
    double worst_margin = 1e9;
    for (int trial = 0; trial < 200000 && certified < 300; ++trial) {
      Genome g = random(rng);
      for (int m = 0; m < 3; ++m) g = mutate(g, rng);
      const Lattice lat = decode(g);
      if (!valid(lat)) continue;
      const FlatInfo f = analyze(lat, kScreenGrid);
      if (f.width > 0.05 || f.gap <= 0.05) continue;
      ++certified;
      for (int k = 0; k < 2000; ++k, ++points) {
        const auto e = energies_at(lat, rng.uniform(0, 2 * M_PI), rng.uniform(0, 2 * M_PI));
        double gap = 1e9;
        if (f.band > 0) gap = std::min(gap, e[f.band] - e[f.band - 1]);
        if (f.band + 1 < lat.n) gap = std::min(gap, e[f.band + 1] - e[f.band]);
        worst_margin = std::min(worst_margin, gap - f.gap);
        if (gap < f.gap - 1e-9) ++violations;
      }
    }
    char buf[256];
    std::snprintf(buf, sizeof buf,
                  "Teste de estresse do certificado: %d redes aleatórias com banda (quase) plana e gap certificado, "
                  "%ld pontos k contínuos sorteados, %d violações (menor folga observada: %.3f).",
                  certified, points, violations, worst_margin);
    return buf;
  }

  static Lattice decode(const Genome& g) {
    Lattice lat;
    lat.n = 1 + g[0] % kMaxSites;
    const double onsite_value[3] = {0, -1, 1};
    for (int s = 0; s < kMaxSites; ++s) lat.onsite[s] = onsite_value[g[kOnsite + s] % 3];
    for (size_t b = 0; b < bonds().size(); ++b) {
      const Bond& bond = bonds()[b];
      const int v = g[kBonds + b] % 3;
      if (v == 0 || bond.i >= lat.n || bond.j >= lat.n) continue;
      lat.hops.push_back({bond, v == 1 ? -1.0 : 1.0});
    }
    return lat;
  }

  // Rede infinita conexa e realmente 2D; todo sítio com ≥ 2 ligações (sem átomos soltos/dímeros triviais).
  static bool valid(const Lattice& lat) {
    std::vector<int> degree(lat.n, 0);
    for (const auto& [b, t] : lat.hops) ++degree[b.i], ++degree[b.j];
    for (int d : degree)
      if (d < 2) return false;
    // BFS no grafo quociente guardando a posição de cada sítio; ciclos geram vetores de rede.
    std::vector<int> px(lat.n, 0), py(lat.n, 0), seen(lat.n, 0);
    std::vector<std::pair<int, int>> cycles;
    seen[0] = 1;
    for (bool changed = true; changed;) {
      changed = false;
      for (const auto& [b, t] : lat.hops) {
        const int u = b.i, v = b.j;
        if (seen[u] && !seen[v]) {
          px[v] = px[u] + b.dx, py[v] = py[u] + b.dy, seen[v] = 1, changed = true;
        } else if (!seen[u] && seen[v]) {
          px[u] = px[v] - b.dx, py[u] = py[v] - b.dy, seen[u] = 1, changed = true;
        }
      }
    }
    for (int s = 0; s < lat.n; ++s)
      if (!seen[s]) return false;
    for (const auto& [b, t] : lat.hops) {
      const int cx = px[b.i] + b.dx - px[b.j], cy = py[b.i] + b.dy - py[b.j];
      if (cx || cy) cycles.push_back({cx, cy});
    }
    for (size_t a = 0; a < cycles.size(); ++a)
      for (size_t c = a + 1; c < cycles.size(); ++c)
        if (cycles[a].first * cycles[c].second - cycles[a].second * cycles[c].first != 0) return true;
    return false;
  }

  static std::vector<std::complex<double>> hamiltonian(const Lattice& lat, double kx, double ky, double shift = 0) {
    const int n = lat.n;
    std::vector<std::complex<double>> h(n * n, 0.0);
    for (int s = 0; s < n; ++s) h[s * n + s] += lat.onsite[s] - shift;
    for (const auto& [b, t] : lat.hops) {
      const auto phase = std::polar(1.0, kx * b.dx + ky * b.dy);
      h[b.i * n + b.j] += t * phase;
      h[b.j * n + b.i] += t * std::conj(phase);
    }
    return h;
  }

  static std::vector<double> energies_at(const Lattice& lat, double kx, double ky) {
    std::complex<double> h[kMaxSites * kMaxSites] = {};
    const int n = lat.n;
    for (int s = 0; s < n; ++s) h[s * n + s] = lat.onsite[s];
    for (const auto& [b, t] : lat.hops) {
      const auto phase = std::polar(1.0, kx * b.dx + ky * b.dy);
      h[b.i * n + b.j] += t * phase;
      h[b.j * n + b.i] += t * std::conj(phase);
    }
    std::vector<double> e(n);
    linalg::small_hermitian_eigenvalues(h, n, e.data());
    return e;
  }

  // Grade grid×grid; de cada par (k, −k) só um é diagonalizado (mesmas energias).
  static Bands compute_bands(const Lattice& lat, int grid) {
    Bands out;
    out.n = lat.n;
    for (int a = 0; a < grid; ++a)
      for (int c = 0; c < grid; ++c) {
        const int na = (grid - a) % grid, nc = (grid - c) % grid;
        if (na * grid + nc < a * grid + c) continue;  // o parceiro −k já foi calculado
        out.e.push_back(energies_at(lat, 2 * M_PI * a / grid, 2 * M_PI * c / grid));
      }
    return out;
  }

  // det(H(k) − E) como polinômio trigonométrico exato, via DFT de (2n+1)² amostras.
  static TrigPoly det_poly(const Lattice& lat, double energy) {
    TrigPoly p;
    p.deg = lat.n;
    const int m = 2 * p.deg + 1;
    std::vector<double> samples(m * m);
    for (int s = 0; s < m; ++s)
      for (int t = 0; t < m; ++t)
        samples[s * m + t] =
            linalg::determinant(hamiltonian(lat, 2 * M_PI * s / m, 2 * M_PI * t / m, energy), lat.n).real();
    p.sample_min = *std::min_element(samples.begin(), samples.end());
    p.sample_max = *std::max_element(samples.begin(), samples.end());
    // DFT separável (primeiro em y, depois em x) com fatores e^{−2πi·j/m} pré-calculados.
    std::vector<std::complex<double>> twiddle(m);
    for (int j = 0; j < m; ++j) twiddle[j] = std::polar(1.0, -2 * M_PI * j / m);
    auto w = [&](int freq, int pos) { return twiddle[(((freq * pos) % m) + m) % m]; };
    std::vector<std::complex<double>> half(m * m);
    for (int s = 0; s < m; ++s)
      for (int b = -p.deg; b <= p.deg; ++b) {
        std::complex<double> acc = 0;
        for (int t = 0; t < m; ++t) acc += samples[s * m + t] * w(b, t);
        half[s * m + b + p.deg] = acc;
      }
    p.c.assign(m * m, 0.0);
    for (int a = -p.deg; a <= p.deg; ++a)
      for (int b = -p.deg; b <= p.deg; ++b) {
        std::complex<double> acc = 0;
        for (int s = 0; s < m; ++s) acc += half[s * m + b + p.deg] * w(a, s);
        const auto coef = acc / static_cast<double>(m * m);
        p.c[(a + p.deg) * m + b + p.deg] = coef;
        p.l1 += std::abs(coef);
        p.cubic += std::abs(coef) * std::pow(a * a + b * b, 1.5);
      }
    return p;
  }

  // Prova que det(H(k) − E) ≠ 0 para todo k (nenhuma banda cruza a energia E).
  // Branch-and-bound com limite de Taylor rigoroso em cada célula (ver abaixo).
  // Retorna false se achar um zero (mudança de sinal ou |f| ≈ 0) ou se não conseguir provar.
  static bool never_zero(const TrigPoly& p) {
    const double tol = 1e-10 * p.l1;
    if (p.l1 < 1e-12) return false;  // identicamente zero
    if (p.sample_min <= tol && p.sample_max >= -tol) return false;  // muda de sinal (ou zera) nas amostras
    struct Cell { double x, y, h; };
    std::vector<Cell> stack;
    stack.reserve(256);
    const int start = 4;
    for (int i = 0; i < start; ++i)
      for (int j = 0; j < start; ++j)
        stack.push_back({2 * M_PI * (i + 0.5) / start, 2 * M_PI * (j + 0.5) / start, M_PI / start});
    int sign = 0;
    long visited = 0;
    while (!stack.empty()) {
      const Cell cell = stack.back();
      stack.pop_back();
      if (++visited > 200000) return false;  // não conseguiu provar: conservador
      const auto l = p.eval(cell.x, cell.y);
      if (std::fabs(l.f) <= tol) return false;
      const int s = l.f > 0 ? 1 : -1;
      if (sign == 0) sign = s;
      if (s != sign) return false;  // mudou de sinal: há um zero no caminho (continuidade)
      // Taylor de 2ª ordem com hessiana local exata + resto cúbico global:
      // |f(c+d)| ≥ |f| − |∇f|·r − ½‖H‖·r² − T₃·r³/6, com ‖H‖ o maior |autovalor| da hessiana 2×2.
      const double r = cell.h * M_SQRT2;
      const double hmean = 0.5 * (l.fxx + l.fyy), hdiff = std::hypot(0.5 * (l.fxx - l.fyy), l.fxy);
      const double hnorm = std::fabs(hmean) + hdiff;
      if (std::fabs(l.f) - std::hypot(l.fx, l.fy) * r - 0.5 * hnorm * r * r - p.cubic * r * r * r / 6 > 0) continue;
      if (cell.h < 1e-9) return false;
      const double q = cell.h / 2;
      for (int dx : {-1, 1})
        for (int dy : {-1, 1}) stack.push_back({cell.x + dx * q, cell.y + dy * q, q});
    }
    return true;
  }

  // Limite inferior RIGOROSO do gap acima de (base) com gap aparente G: prova que nenhuma banda
  // passa por dois níveis base + φ1·G e base + φ2·G; então o gap real é ≥ (φ2 − φ1)·G.
  // Se o lado de baixo da janela é uma banda EXATAMENTE plana (energia conhecida), um nível basta:
  // nada cruza base + φ·G  ⇒  gap ≥ φ·G.
  static double certify_side(const Lattice& lat, double base, double apparent, bool exact_base) {
    if (apparent < 1e-3) return 0;
    if (exact_base) {
      for (double phi : {0.8, 0.5, 0.2})
        if (never_zero(det_poly(lat, base + phi * apparent))) return phi * apparent;
      return 0;
    }
    const double levels[3][2] = {{0.1, 0.9}, {0.3, 0.7}, {0.45, 0.55}};
    for (const auto& l : levels)
      if (never_zero(det_poly(lat, base + l[0] * apparent)) && never_zero(det_poly(lat, base + l[1] * apparent)))
        return (l[1] - l[0]) * apparent;
    return 0;
  }

  static double certify_side_down(const Lattice& lat, double flat_energy, double apparent) {
    if (apparent < 1e-3) return 0;
    for (double phi : {0.8, 0.5, 0.2})
      if (never_zero(det_poly(lat, flat_energy - phi * apparent))) return phi * apparent;
    return 0;
  }

  static double band_min(const Bands& b, int band) {
    double v = 1e9;
    for (const auto& e : b.e) v = std::min(v, e[band]);
    return v;
  }
  static double band_max(const Bands& b, int band) {
    double v = -1e9;
    for (const auto& e : b.e) v = std::max(v, e[band]);
    return v;
  }

  // Banda mais plana (menor largura na grade) e seu gap aparente (separação global de energia).
  static FlatInfo flattest(const Bands& b) {
    FlatInfo best;
    for (int band = 0; band < b.n; ++band) {
      const double lo = band_min(b, band), hi = band_max(b, band);
      if (hi - lo >= best.width) continue;
      double gap = b.n > 1 ? 1e9 : 0;
      if (band > 0) gap = std::min(gap, lo - band_max(b, band - 1));
      if (band + 1 < b.n) gap = std::min(gap, band_min(b, band + 1) - hi);
      best.band = band, best.width = hi - lo, best.energy = 0.5 * (lo + hi), best.grid_gap = std::max(0.0, gap);
    }
    return best;
  }

  // Avaliação completa. Dispersivas: só a largura conta (score ≤ ~1.3).
  // Planas/quase planas: largura na grade fina + gap certificado.
  static FlatInfo analyze(const Lattice& lat, int screen_grid, int fine_grid = kGrid) {
    FlatInfo f = flattest(compute_bands(lat, screen_grid));
    if (f.width <= 0.05) {
      const Bands bands = compute_bands(lat, fine_grid);
      f = flattest(bands);
      if (f.width <= 0.05) {
        // Janelas acima e abaixo da banda. Para banda exatamente plana, a janela é medida a partir
        // da própria banda (energia conhecida): abaixo, espelhamos com −H para reusar o mesmo teste.
        const bool exact = f.width < 1e-6;
        double gap = lat.n > 1 ? 1e9 : 0;
        if (f.band > 0) {
          const double below = band_max(bands, f.band - 1), bottom = band_min(bands, f.band);
          gap = std::min(gap, exact ? certify_side_down(lat, bottom, bottom - below)
                                    : certify_side(lat, below, bottom - below, false));
        }
        if (f.band + 1 < lat.n) {
          const double top = band_max(bands, f.band);
          gap = std::min(gap, certify_side(lat, top, band_min(bands, f.band + 1) - top, exact));
        }
        f.gap = gap;
        f.score = -std::log10(f.width + 1e-6) + 3.0 * std::min(f.gap, 1.0);
        return f;
      }
    }
    f.gap = 0;
    f.score = -std::log10(f.width + 1e-6);
    return f;
  }

  static double min_energy(const Bands& b) {
    double v = 1e9;
    for (const auto& e : b.e) v = std::min(v, e.front());
    return v;
  }
  static double max_energy(const Bands& b) {
    double v = -1e9;
    for (const auto& e : b.e) v = std::max(v, e.back());
    return v;
  }

  static std::string classify(const Lattice& lat, const FlatInfo& f) {
    if (f.width > 0.05) return "dispersiva";
    std::string kind = f.width < 1e-6 ? "banda plana" : "quase plana";
    // "isolada" = gap PROVADO ≥ 0.05; "nao isolada" = não foi possível provar (inclui as que se tocam).
    kind += f.gap > 0.05 ? " isolada" : " nao isolada";
    if (f.width < 1e-6) {
      // Origem: desbalanço de subredes (bipartida, sem energias locais, E=0) ou interferência.
      bool bipartite = std::abs(f.energy) < 1e-6;
      for (int s = 0; s < lat.n && bipartite; ++s) bipartite = lat.onsite[s] == 0;
      std::vector<int> color(lat.n, -1);
      if (bipartite) {
        color[0] = 0;
        for (int pass = 0; pass < lat.n + 1; ++pass)
          for (const auto& [b, t] : lat.hops) {
            if (color[b.i] >= 0 && color[b.j] < 0) color[b.j] = 1 - color[b.i];
            if (color[b.j] >= 0 && color[b.i] < 0) color[b.i] = 1 - color[b.j];
          }
        for (const auto& [b, t] : lat.hops)
          if (color[b.i] == color[b.j]) bipartite = false;
      }
      int na = 0;
      for (int c : color) na += (c == 0);
      if (bipartite && 2 * na != lat.n)
        kind += " [desbalanco de subredes]";
      else
        kind += " [interferencia]";
    }
    char buf[32];
    std::snprintf(buf, sizeof buf, ", %d sitios", lat.n);
    return kind + buf;
  }
};

}  // namespace

std::unique_ptr<Lab> make_materials_lab() { return std::make_unique<MaterialsLab>(); }
