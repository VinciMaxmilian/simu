// Fase 2 — Química: hidrocarbonetos policíclicos benzenoides (grafenos moleculares).
//
// Física: teoria de Hückel dos elétrons π (tight-binding no grafo de carbonos) — a mecânica
// quântica de elétrons + núcleos que o Modelo Padrão fornece para a química, na aproximação
// mais simples que ainda é quantitativamente útil para essa família de moléculas.
//
// Objetivo de pesquisa: "carbono magnético". Modos de energia zero (orbitais não ligantes)
// abrigam elétrons desemparelhados → moléculas com spin/magnetismo feitas só de C e H.
//   η          = nº de modos zero = N − 2·posto(B), B = matriz de biadjacência das subredes
//   desbalanço = |N_A − N_B| (Teorema de Lieb/Longuet-Higgins: garante η ≥ desbalanço)
//   oculto     = η − desbalanço: modos zero que o desbalanço NÃO explica (ex.: cálice de Clar, 2021)
// Os "ocultos" são os mais interessantes: magnetismo por interferência topológica.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <queue>

#include "core/linalg.hpp"
#include "labs/labs.hpp"

namespace {

constexpr int R = 4;          // raio da grade hexagonal disponível
constexpr int kMaxHex = 16;   // tamanho máximo da molécula (anéis)

struct Cell { int q, r; };

// Vértices de um hexágono "pointy-top" em coordenadas inteiras (x em √3/2, y em 1/2).
constexpr int kCorner[6][2] = {{0, 2}, {1, 1}, {1, -1}, {0, -2}, {-1, -1}, {-1, 1}};
constexpr int kNeighbor[6][2] = {{1, 0}, {0, 1}, {-1, 1}, {-1, 0}, {0, -1}, {1, -1}};

struct Molecule {
  int atoms = 0, na = 0, nb = 0, hydrogens = 0, hexes = 0;
  std::vector<std::pair<int, int>> bonds;
  std::vector<int> sublattice;  // 0 = A, 1 = B
  std::vector<std::pair<int, int>> position;  // (x, y) inteiros: x em unidades de √3/2·d, y em d/2
};

class ChemistryLab : public Lab {
 public:
  ChemistryLab() {
    for (int q = -R; q <= R; ++q)
      for (int r = -R; r <= R; ++r)
        if (std::abs(q + r) <= R) {
          index_[{q, r}] = static_cast<int>(cells_.size());
          cells_.push_back({q, r});
        }
  }

  std::string id() const override { return "quimica"; }
  std::string title() const override { return "Química: benzenoides e magnetismo de carbono (Hückel)"; }
  std::string objective() const override {
    return "maximizar modos de energia zero (elétrons desemparelhados), com peso dobrado para os 'ocultos' "
           "(não explicados pelo desbalanço de subredes); até 16 anéis";
  }

  Genome random(Rng& rng) const override {
    Genome g(cells_.size(), 0);
    const int size = 1 + rng.below(kMaxHex);
    g[index_.at({0, 0})] = 1;
    for (int n = 1; n < size; ++n) {
      auto frontier = frontier_of(g);
      if (frontier.empty()) break;
      g[frontier[rng.below(static_cast<int>(frontier.size()))]] = 1;
    }
    return g;
  }

  Genome mutate(const Genome& g, Rng& rng) const override {
    Genome m = g;
    const int ops = 1 + rng.below(2);
    for (int o = 0; o < ops; ++o) {
      const int size = count(m);
      const bool add = size < 2 || (size < kMaxHex && rng.chance(0.5));
      if (add) {
        auto frontier = frontier_of(m);
        if (!frontier.empty()) m[frontier[rng.below(static_cast<int>(frontier.size()))]] = 1;
      } else {
        // Remove um anel sem desconectar a molécula.
        std::vector<int> occupied;
        for (int i = 0; i < static_cast<int>(m.size()); ++i)
          if (m[i]) occupied.push_back(i);
        for (int attempt = 0; attempt < 8; ++attempt) {
          const int i = occupied[rng.below(static_cast<int>(occupied.size()))];
          m[i] = 0;
          if (connected(m)) break;
          m[i] = 1;
        }
      }
    }
    return m;
  }

  Genome crossover(const Genome& a, const Genome& b, Rng& rng) const override {
    // Corta o plano por uma reta aleatória: metade de cada "pai", fica o maior pedaço conexo.
    Genome c(cells_.size(), 0);
    const double angle = rng.uniform(0, 2 * M_PI);
    const double dx = std::cos(angle), dy = std::sin(angle);
    for (size_t i = 0; i < cells_.size(); ++i) {
      const double x = std::sqrt(3.0) * (cells_[i].q + cells_[i].r / 2.0), y = 1.5 * cells_[i].r;
      c[i] = (x * dx + y * dy >= 0) ? a[i] : b[i];
    }
    c = largest_component(c);
    Rng local(rng.next());
    while (count(c) > kMaxHex) c = mutate_remove(c, local);
    if (count(c) == 0) c = a;
    return c;
  }

  Evaluation evaluate(const Genome& g) const override {
    Evaluation e;
    const int hexes = count(g);
    if (hexes == 0 || hexes > kMaxHex || !connected(g)) return e;
    const Molecule mol = build(g);
    const int eta = zero_modes(mol);
    const int imbalance = std::abs(mol.na - mol.nb);
    const int hidden = eta - imbalance;
    e.valid = true;
    e.score = eta + 2.0 * hidden - 0.02 * hexes;
    char buf[64];
    if (eta == 0)
      std::snprintf(buf, sizeof buf, "camada fechada (eta=0)");
    else
      std::snprintf(buf, sizeof buf, "eta=%d (desbalanco=%d, oculto=%d)", eta, imbalance, hidden);
    e.label = buf;
    return e;
  }

  // A mesma molécula transladada na grade é a mesma física: chave = células relativas à menor (q, r).
  std::string canonical(const Genome& g) const override {
    int mq = 1 << 30, mr = 1 << 30;
    for (size_t i = 0; i < g.size(); ++i)
      if (g[i]) mr = std::min(mr, cells_[i].r);
    for (size_t i = 0; i < g.size(); ++i)
      if (g[i] && cells_[i].r == mr) mq = std::min(mq, cells_[i].q);
    std::string key;
    for (size_t i = 0; i < g.size(); ++i)
      if (g[i]) key += {static_cast<char>(cells_[i].q - mq + 16), static_cast<char>(cells_[i].r - mr + 16)};
    return key;
  }

  std::string describe(const Genome& g) const override {
    const Molecule mol = build(g);
    const int eta = zero_modes(mol);
    char buf[256];
    std::snprintf(buf, sizeof buf, "C%dH%d, %d anéis, subredes %d/%d, modos zero η=%d, gap HOMO-LUMO=%.3f|β|",
                  mol.atoms, mol.hydrogens, mol.hexes, mol.na, mol.nb, eta, gap(mol));
    return std::string(buf) + "\n" + picture(g);
  }

  // Tudo para reconstruir a molécula: anéis (coordenadas axiais), átomos (Å, ligação C–C = 1.42 Å),
  // ligações, subredes, e o que o Hückel calculou.
  std::string export_json(const Genome& g) const override {
    const Molecule mol = build(g);
    std::string s = "{\"anel_coordenadas_axiais_q_r\": [";
    bool first = true;
    for (size_t i = 0; i < g.size(); ++i)
      if (g[i]) {
        s += (first ? "" : ", ") + std::string("[") + std::to_string(cells_[i].q) + ", " + std::to_string(cells_[i].r) + "]";
        first = false;
      }
    char buf[160];
    std::snprintf(buf, sizeof buf, "], \"formula\": \"C%dH%d\", \"aneis\": %d, \"eta\": %d, \"n_A\": %d, \"n_B\": %d, "
                  "\"gap_huckel_beta\": %.6f, ", mol.atoms, mol.hydrogens, mol.hexes, zero_modes(mol), mol.na, mol.nb, gap(mol));
    s += buf;
    s += "\"carbonos\": [";
    constexpr double d = 1.42;
    for (int i = 0; i < mol.atoms; ++i) {
      std::snprintf(buf, sizeof buf, "%s{\"id\": %d, \"x_A\": %.5f, \"y_A\": %.5f, \"subrede\": \"%c\"}", i ? ", " : "", i,
                    mol.position[i].first * std::sqrt(3.0) / 2 * d, mol.position[i].second * 0.5 * d,
                    mol.sublattice[i] == 0 ? 'A' : 'B');
      s += buf;
    }
    s += "], \"ligacoes\": [";
    for (size_t k = 0; k < mol.bonds.size(); ++k)
      s += (k ? ", [" : "[") + std::to_string(mol.bonds[k].first) + ", " + std::to_string(mol.bonds[k].second) + "]";
    s += "]}";
    return s;
  }

  std::vector<std::string> validate() const override {
    struct Case { const char* name; std::vector<Cell> cells; const char* formula; double gap; int eta; };
    const std::vector<Case> cases = {
        {"benzeno", {{0, 0}}, "C6H6", 2.000, 0},
        {"naftaleno", {{0, 0}, {1, 0}}, "C10H8", 1.236, 0},
        {"antraceno", {{0, 0}, {1, 0}, {2, 0}}, "C14H10", 0.828, 0},
        {"fenantreno", {{0, 0}, {1, 0}, {1, 1}}, "C14H10", 1.210, 0},
        {"fenalenil (radical)", {{0, 0}, {1, 0}, {0, 1}}, "C13H9", 0.000, 1},
        {"[3]triangulena (sintetizada em 2017)", {{0, 0}, {1, 0}, {2, 0}, {0, 1}, {1, 1}, {0, 2}}, "C22H12", 0.000, 2},
    };
    std::vector<std::string> out;
    out.push_back("| molécula | fórmula (real → modelo) | gap Hückel (literatura → modelo) | η (esperado → modelo) |");
    out.push_back("|---|---|---|---|");
    for (const auto& c : cases) {
      Genome g(cells_.size(), 0);
      for (const auto& cell : c.cells) g[index_.at({cell.q, cell.r})] = 1;
      const Molecule mol = build(g);
      char buf[256];
      std::snprintf(buf, sizeof buf, "| %s | %s → C%dH%d | %.3f → %.3f | %d → %d |", c.name, c.formula, mol.atoms,
                    mol.hydrogens, c.gap, gap(mol), c.eta, zero_modes(mol));
      out.push_back(buf);
    }
    return out;
  }

 private:
  static int count(const Genome& g) {
    int n = 0;
    for (auto b : g) n += b;
    return n;
  }

  std::vector<int> neighbors(int i) const {
    std::vector<int> out;
    for (const auto& d : kNeighbor) {
      auto it = index_.find({cells_[i].q + d[0], cells_[i].r + d[1]});
      if (it != index_.end()) out.push_back(it->second);
    }
    return out;
  }

  std::vector<int> frontier_of(const Genome& g) const {
    std::vector<int> out;
    for (int i = 0; i < static_cast<int>(g.size()); ++i) {
      if (g[i]) continue;
      for (int n : neighbors(i))
        if (g[n]) { out.push_back(i); break; }
    }
    return out;
  }

  std::vector<int> component(const Genome& g, int start, std::vector<uint8_t>& seen) const {
    std::vector<int> comp;
    std::queue<int> todo;
    todo.push(start);
    seen[start] = 1;
    while (!todo.empty()) {
      const int i = todo.front();
      todo.pop();
      comp.push_back(i);
      for (int n : neighbors(i))
        if (g[n] && !seen[n]) seen[n] = 1, todo.push(n);
    }
    return comp;
  }

  bool connected(const Genome& g) const {
    std::vector<uint8_t> seen(g.size(), 0);
    for (int i = 0; i < static_cast<int>(g.size()); ++i)
      if (g[i]) return static_cast<int>(component(g, i, seen).size()) == count(g);
    return false;
  }

  Genome largest_component(const Genome& g) const {
    std::vector<uint8_t> seen(g.size(), 0);
    std::vector<int> best;
    for (int i = 0; i < static_cast<int>(g.size()); ++i)
      if (g[i] && !seen[i]) {
        auto comp = component(g, i, seen);
        if (comp.size() > best.size()) best = comp;
      }
    Genome out(g.size(), 0);
    for (int i : best) out[i] = 1;
    return out;
  }

  Genome mutate_remove(const Genome& g, Rng& rng) const {
    Genome m = g;
    std::vector<int> occupied;
    for (int i = 0; i < static_cast<int>(m.size()); ++i)
      if (m[i]) occupied.push_back(i);
    for (int attempt = 0; attempt < 32; ++attempt) {
      const int i = occupied[rng.below(static_cast<int>(occupied.size()))];
      m[i] = 0;
      if (connected(m)) return m;
      m[i] = 1;
    }
    m[occupied[rng.below(static_cast<int>(occupied.size()))]] = 0;
    return largest_component(m);
  }

  Molecule build(const Genome& g) const {
    Molecule mol;
    std::map<std::pair<int, int>, int> atom;
    std::map<std::pair<int, int>, int> bond_set;
    auto atom_id = [&](int x, int y) {
      auto [it, inserted] = atom.emplace(std::pair{x, y}, mol.atoms);
      if (inserted) {
        ++mol.atoms;
        mol.position.push_back({x, y});
        mol.sublattice.push_back((((y % 3) + 3) % 3) == 2 ? 0 : 1);
      }
      return it->second;
    };
    for (int i = 0; i < static_cast<int>(g.size()); ++i) {
      if (!g[i]) continue;
      ++mol.hexes;
      const int cx = 2 * cells_[i].q + cells_[i].r, cy = 3 * cells_[i].r;
      int ids[6];
      for (int k = 0; k < 6; ++k) ids[k] = atom_id(cx + kCorner[k][0], cy + kCorner[k][1]);
      for (int k = 0; k < 6; ++k) {
        int a = ids[k], b = ids[(k + 1) % 6];
        if (a > b) std::swap(a, b);
        bond_set[{a, b}] = 1;
      }
    }
    for (const auto& [ab, _] : bond_set) mol.bonds.push_back(ab);
    std::vector<int> degree(mol.atoms, 0);
    for (auto [a, b] : mol.bonds) ++degree[a], ++degree[b];
    for (int d : degree) mol.hydrogens += (d == 2);
    for (int s : mol.sublattice) (s == 0 ? mol.na : mol.nb)++;
    return mol;
  }

  static int zero_modes(const Molecule& mol) {
    std::vector<int> pos(mol.atoms);
    int ia = 0, ib = 0;
    for (int i = 0; i < mol.atoms; ++i) pos[i] = mol.sublattice[i] == 0 ? ia++ : ib++;
    std::vector<int64_t> b(static_cast<size_t>(mol.na) * mol.nb, 0);
    for (auto [u, v] : mol.bonds) {
      if (mol.sublattice[u] == 1) std::swap(u, v);
      b[pos[u] * mol.nb + pos[v]] = 1;
    }
    return mol.atoms - 2 * linalg::rank_exact(b, mol.na, mol.nb);
  }

  static double gap(const Molecule& mol) {
    const int n = mol.atoms;
    std::vector<double> h(n * n, 0.0);
    for (auto [a, b] : mol.bonds) h[a * n + b] = h[b * n + a] = 1.0;  // unidades de |β|
    auto eig = linalg::symmetric_eigenvalues(h, n);
    // Energia E = α + xβ com β < 0: x maior = mais ligante. n elétrons π ocupam n/2 orbitais.
    std::sort(eig.rbegin(), eig.rend());
    const int homo = (n + 1) / 2 - 1;
    if (n % 2 == 1) return 0.0;  // radical com número ímpar de elétrons: SOMO semipreenchido
    return eig[homo] - eig[homo + 1];
  }

  std::string picture(const Genome& g) const {
    std::string out;
    for (int r = -R; r <= R; ++r) {
      std::string row(std::abs(r), ' ');
      bool any = false;
      for (int q = -R; q <= R; ++q) {
        auto it = index_.find({q, r});
        if (it == index_.end()) continue;
        row += g[it->second] ? "⬢ " : "· ";
        any = any || g[it->second];
      }
      if (any) out += "    " + row + "\n";
    }
    return out;
  }

  std::vector<Cell> cells_;
  std::map<std::pair<int, int>, int> index_;
};

}  // namespace

std::unique_ptr<Lab> make_chemistry_lab() { return std::make_unique<ChemistryLab>(); }
