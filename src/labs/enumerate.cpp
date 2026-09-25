// Enumeração exaustiva de benzenoides (poli-hexágonos) até h anéis: a prova que a busca heurística
// das civilizações não dá. Algoritmo de Redelmeier na rede hexagonal (formas fixas); cada forma livre
// é avaliada uma vez, no representante canônico entre as 12 simetrias do hexágono.
// Para cada forma: fórmula, subredes, modos zero (posto exato) e planaridade (H–H ≥ 1.5 Å).
#include "labs/enumerate.hpp"

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

#include "core/linalg.hpp"

namespace {

constexpr int kNeighbor[6][2] = {{1, 0}, {0, 1}, {-1, 1}, {-1, 0}, {0, -1}, {1, -1}};
constexpr int kCorner[6][2] = {{0, 2}, {1, 1}, {1, -1}, {0, -2}, {-1, -1}, {-1, 1}};
constexpr int kOff = 32, kW = 64;  // grade de trabalho: células (q, r) com |q|, |r| < 32

using Cells = std::vector<std::pair<int, int>>;

// 12 simetrias do hexágono em coordenadas cúbicas (x=q, z=r, y=−q−r): rotações e reflexões.
std::pair<int, int> transform(int q, int r, int k) {
  int x = q, z = r, y = -q - r;
  for (int i = 0; i < k % 6; ++i) {  // rotação de 60°: (x, y, z) → (−z, −x, −y)
    const int nx = -z, ny = -x, nz = -y;
    x = nx, y = ny, z = nz;
  }
  if (k >= 6) std::swap(x, z);  // reflexão
  return {x, z};
}

Cells normalize(Cells c) {
  std::sort(c.begin(), c.end(), [](auto a, auto b) { return a.second != b.second ? a.second < b.second : a.first < b.first; });
  const int q0 = c[0].first, r0 = c[0].second;
  for (auto& [q, r] : c) q -= q0, r -= r0;
  return c;
}

bool is_canonical(const Cells& cells) {
  const Cells base = normalize(cells);
  for (int k = 1; k < 12; ++k) {
    Cells t;
    t.reserve(cells.size());
    for (auto [q, r] : cells) t.push_back(transform(q, r, k));
    if (normalize(t) < base) return false;
  }
  return true;
}

// Avalia uma forma: átomos, subredes, modos zero e menor H–H (mesmas convenções de chemistry.cpp).
BenzenoidInfo evaluate(const Cells& cells) {
  std::map<std::pair<int, int>, int> atom;
  std::vector<std::pair<int, int>> pos;
  std::vector<std::pair<int, int>> bonds;
  for (auto [q, r] : cells) {
    const int cx = 2 * q + r, cy = 3 * r;
    int ids[6];
    for (int k = 0; k < 6; ++k) {
      auto key = std::pair{cx + kCorner[k][0], cy + kCorner[k][1]};
      auto [it, ins] = atom.emplace(key, static_cast<int>(pos.size()));
      if (ins) pos.push_back(key);
      ids[k] = it->second;
    }
    for (int k = 0; k < 6; ++k) bonds.push_back(std::minmax(ids[k], ids[(k + 1) % 6]));
  }
  std::sort(bonds.begin(), bonds.end());
  bonds.erase(std::unique(bonds.begin(), bonds.end()), bonds.end());
  const int n = static_cast<int>(pos.size());
  std::vector<std::vector<int>> nb(n);
  for (auto [a, b] : bonds) nb[a].push_back(b), nb[b].push_back(a);

  BenzenoidInfo info;
  info.rings = static_cast<int>(cells.size());
  info.carbons = n;
  std::vector<int> sub(n), idx(n);
  for (int i = 0; i < n; ++i) {
    sub[i] = (((pos[i].second % 3) + 3) % 3) == 2 ? 0 : 1;
    idx[i] = sub[i] == 0 ? info.n_a++ : info.n_b++;
    info.hydrogens += nb[i].size() == 2;
  }
  // Planaridade primeiro (barata): H–H mínimo na geometria plana.
  constexpr double d = 1.42, ch = 1.09;
  std::vector<std::pair<double, double>> h;
  for (int i = 0; i < n; ++i) {
    if (nb[i].size() != 2) continue;
    const double x = pos[i].first * std::sqrt(3.0) / 2 * d, y = pos[i].second * 0.5 * d;
    double vx = x, vy = y;
    for (int j : nb[i]) vx -= 0.5 * pos[j].first * std::sqrt(3.0) / 2 * d, vy -= 0.5 * pos[j].second * 0.5 * d;
    const double norm = std::hypot(vx, vy);
    h.push_back({x + ch * vx / norm, y + ch * vy / norm});
  }
  info.min_hh = 1e9;
  for (size_t a = 0; a < h.size(); ++a)
    for (size_t b = a + 1; b < h.size(); ++b)
      info.min_hh = std::min(info.min_hh, std::hypot(h[a].first - h[b].first, h[a].second - h[b].second));
  std::vector<int64_t> m(static_cast<size_t>(info.n_a) * info.n_b, 0);
  for (auto [a, b] : bonds) {
    const int u = sub[a] == 0 ? a : b, v = sub[a] == 0 ? b : a;
    m[idx[u] * info.n_b + idx[v]] = 1;
  }
  info.eta = n - 2 * linalg::rank_exact(m, info.n_a, info.n_b);
  info.cells = normalize(cells);
  return info;
}

struct Redelmeier {
  int max_h, split_depth, thread, threads;
  const std::function<void(const BenzenoidInfo&)>* visit;
  std::vector<long>* fixed_counts;
  std::vector<long>* free_counts;
  std::array<uint8_t, kW * kW> used{};  // já considerada (na forma ou na fronteira)
  Cells current;
  long branch = 0;

  static bool allowed(int q, int r) { return r > 0 || (r == 0 && q >= 0); }  // célula inicial é a menor

  void run() {
    Cells untried = {{0, 0}};
    used[kOff * kW + kOff] = 1;
    recurse(untried, 0);
  }

  void recurse(Cells untried, int depth) {
    while (!untried.empty()) {
      const auto cell = untried.back();
      untried.pop_back();
      // Divide o trabalho entre threads na profundidade split_depth.
      if (depth == split_depth && (branch++ % threads) != thread) continue;
      current.push_back(cell);
      const int h = static_cast<int>(current.size());
      if (depth >= split_depth) {
        (*fixed_counts)[h]++;
        if (is_canonical(current)) {
          (*free_counts)[h]++;
          (*visit)(evaluate(current));
        }
      } else if (thread == 0) {
        (*fixed_counts)[h]++;
        if (is_canonical(current)) {
          (*free_counts)[h]++;
          (*visit)(evaluate(current));
        }
      }
      if (h < max_h) {
        Cells next = untried, added;
        for (const auto& dlt : kNeighbor) {
          const int q = cell.first + dlt[0], r = cell.second + dlt[1];
          if (!allowed(q, r) || used[(r + kOff) * kW + q + kOff]) continue;
          used[(r + kOff) * kW + q + kOff] = 1;
          added.push_back({q, r});
          next.push_back({q, r});
        }
        recurse(next, depth + 1);
        for (auto [q, r] : added) used[(r + kOff) * kW + q + kOff] = 0;
      }
      current.pop_back();
    }
  }
};

}  // namespace

EnumerationResult enumerate_benzenoids(int max_h, const std::function<void(const BenzenoidInfo&)>& visit) {
  EnumerationResult out;
  out.fixed.assign(max_h + 1, 0);
  out.free.assign(max_h + 1, 0);
  const int threads = omp_get_max_threads();
  const int split = std::min(4, max_h - 1);
#pragma omp parallel num_threads(threads)
  {
    std::vector<long> fixed(max_h + 1, 0), fre(max_h + 1, 0);
    Redelmeier r{max_h, split, omp_get_thread_num(), threads, nullptr, &fixed, &fre, {}, {}, 0};
    std::function<void(const BenzenoidInfo&)> local = [&](const BenzenoidInfo& b) {
#pragma omp critical(enum_visit)
      visit(b);
    };
    r.visit = &local;
    r.run();
#pragma omp critical(enum_merge)
    for (int h = 0; h <= max_h; ++h) out.fixed[h] += fixed[h], out.free[h] += fre[h];
  }
  return out;
}
