// Fase 4 — Mundo-brinquedo.
//
// Um universo artificial cujas "leis" nós conhecemos por completo, para validar a engrenagem
// civilização ↔ laboratório: sabemos o ótimo global, então dá para medir o quão perto as
// civilizações chegam e quais leis escondidas elas descobrem.
//
// Física do brinquedo: paisagem NK (Kauffman) com L=28 variáveis e epistasia K=5, mais 4 "leis"
// plantadas (padrões de 5 variáveis que dão bônus). As leis podem conflitar entre si.
#include <omp.h>

#include <array>
#include <cmath>
#include <cstdio>

#include "labs/labs.hpp"

namespace {

constexpr int L = 28;
constexpr int K = 5;
constexpr int kLawSize = 5;

struct Law {
  std::string name;
  std::array<int, kLawSize> pos;
  std::array<uint8_t, kLawSize> val;
  double bonus;
};

class ToyLab : public Lab {
 public:
  explicit ToyLab(uint64_t seed) {
    Rng rng(seed);
    for (int i = 0; i < L; ++i) {
      std::array<int, K> nb;
      for (int k = 0; k < K; ++k) {
        int j;
        do { j = rng.below(L); } while (j == i);
        nb[k] = j;
      }
      neighbors_[i] = nb;
      for (auto& v : table_[i]) v = rng.uniform();
    }
    const char* names[] = {"Lei-A", "Lei-B", "Lei-C", "Lei-D"};
    for (int l = 0; l < 4; ++l) {
      Law law{names[l], {}, {}, 0.12};
      for (int k = 0; k < kLawSize; ++k) {
        bool repeated;
        do {  // posições distintas: uma lei nunca pode contradizer a si mesma
          law.pos[k] = rng.below(L);
          repeated = false;
          for (int j = 0; j < k; ++j) repeated = repeated || law.pos[j] == law.pos[k];
        } while (repeated);
        law.val[k] = static_cast<uint8_t>(rng.below(2));
      }
      laws_.push_back(law);
    }
    // Força bruta: 2^28 configurações (em paralelo), para conhecer a verdade.
    optimum_ = -1;
#pragma omp parallel
    {
      double local_best = -1;
      uint32_t local_x = 0;
      Genome g(L);
#pragma omp for schedule(static)
      for (int64_t x = 0; x < (int64_t{1} << L); ++x) {
        for (int i = 0; i < L; ++i) g[i] = (x >> i) & 1;
        const double f = raw(g);
        if (f > local_best) local_best = f, local_x = static_cast<uint32_t>(x);
      }
#pragma omp critical
      if (local_best > optimum_) {
        optimum_ = local_best;
        best_.assign(L, 0);
        for (int i = 0; i < L; ++i) best_[i] = (local_x >> i) & 1;
      }
    }
  }

  std::string id() const override { return "brinquedo"; }
  std::string title() const override { return "Mundo-brinquedo (paisagem NK com leis escondidas)"; }
  std::string objective() const override {
    return "maximizar o desempenho de uma 'tecnologia' de 28 componentes; score = fração do ótimo global";
  }

  Genome random(Rng& rng) const override {
    Genome g(L);
    for (auto& b : g) b = static_cast<uint8_t>(rng.below(2));
    return g;
  }
  Genome mutate(const Genome& g, Rng& rng) const override {
    Genome m = g;
    const int flips = 1 + rng.below(2);
    for (int f = 0; f < flips; ++f) m[rng.below(L)] ^= 1;
    return m;
  }
  Genome crossover(const Genome& a, const Genome& b, Rng& rng) const override {
    Genome c(L);
    for (int i = 0; i < L; ++i) c[i] = rng.chance(0.5) ? a[i] : b[i];
    return c;
  }

  Evaluation evaluate(const Genome& g) const override {
    Evaluation e;
    e.valid = true;
    e.score = raw(g) / optimum_;
    std::string label;
    for (const auto& law : laws_)
      if (satisfies(g, law)) label += (label.empty() ? "" : "+") + law.name;
    e.label = label.empty() ? "nenhuma lei" : label;
    return e;
  }

  std::string describe(const Genome& g) const override {
    std::string bits;
    for (auto b : g) bits += static_cast<char>('0' + b);
    char buf[128];
    std::snprintf(buf, sizeof buf, "config %s — %.2f%% do ótimo", bits.c_str(), 100.0 * raw(g) / optimum_);
    return buf;
  }

  std::vector<std::string> validate() const override {
    std::vector<std::string> out;
    out.push_back("Ótimo global (força bruta sobre 2^28 configurações): " + describe(best_));
    std::string opt_label = evaluate(best_).label;
    out.push_back("Leis satisfeitas pelo ótimo global: " + opt_label);
    for (const auto& law : laws_) {
      std::string p;
      for (int k = 0; k < kLawSize; ++k) p += "x" + std::to_string(law.pos[k]) + "=" + std::to_string(law.val[k]) + " ";
      out.push_back(law.name + ": " + p + "(bônus " + std::to_string(law.bonus).substr(0, 4) + ")");
    }
    return out;
  }

 private:
  static bool satisfies(const Genome& g, const Law& law) {
    for (int k = 0; k < kLawSize; ++k)
      if (g[law.pos[k]] != law.val[k]) return false;
    return true;
  }
  double raw(const Genome& g) const {
    double f = 0;
    for (int i = 0; i < L; ++i) {
      int idx = g[i];
      for (int k = 0; k < K; ++k) idx = (idx << 1) | g[neighbors_[i][k]];
      f += table_[i][idx];
    }
    f /= L;
    for (const auto& law : laws_)
      if (satisfies(g, law)) f += law.bonus;
    return f;
  }

  std::array<std::array<int, K>, L> neighbors_;
  std::array<std::array<double, 1 << (K + 1)>, L> table_;
  std::vector<Law> laws_;
  double optimum_;
  Genome best_;
};

}  // namespace

std::unique_ptr<Lab> make_toy_lab(uint64_t seed) { return std::make_unique<ToyLab>(seed); }
