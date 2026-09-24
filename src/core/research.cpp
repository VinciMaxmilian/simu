#include "core/research.hpp"

#include <algorithm>
#include <cmath>
#include <set>
#include <unordered_map>

std::vector<Culture> default_cultures() {
  return {
      {"Empiristas", 0.50, 0.10, 0.15},        // testam de tudo
      {"Incrementalistas", 0.03, 0.05, 0.15},  // refinam o que já funciona
      {"Evolucionistas", 0.10, 0.45, 0.15},    // recombinam ideias
      {"Isolados", 0.15, 0.20, 0.0},           // equilibrados, mas sem contato
  };
}

namespace {

struct Entry {
  Genome g;
  Evaluation e;
};

struct Civ {
  const Culture* culture;
  std::vector<Entry> archive;              // memória de trabalho
  std::map<std::string, Entry> catalog;    // melhor exemplo de cada classe conhecida
  Entry best;
  bool has_best = false;
  int dark_until = -1;
  long experiments = 0;
  int catastrophes = 0;
};

class Engine {
 public:
  Engine(const Lab& lab, const std::vector<Culture>& cultures, const ResearchConfig& cfg, uint64_t seed,
         int run_id)
      : lab_(lab), cfg_(cfg), rng_(seed), run_id_(run_id) {
    for (const auto& c : cultures) civs_.push_back(Civ{&c, {}, {}, {}, false, -1, 0, 0});
  }

  RunResult run() {
    const int n = static_cast<int>(civs_.size());
    for (int year = 0; year <= cfg_.years; ++year) {
      for (int c = 0; c < n; ++c) research_year(c, year);
      for (int c = 0; c < n; ++c) exchange(c, year);
      for (int c = 0; c < n; ++c) maybe_catastrophe(c, year);
      if (year % cfg_.log_every == 0 || year == cfg_.years)
        for (int c = 0; c < n; ++c) {
          const Civ& civ = civs_[c];
          out_.timeline.push_back({run_id_, year, c, civ.has_best ? civ.best.e.score : NAN,
                                   static_cast<int>(civ.catalog.size()), civ.experiments});
        }
    }
    for (const Civ& civ : civs_) {
      out_.best_genome.push_back(civ.best.g);
      out_.best_score.push_back(civ.has_best ? civ.best.e.score : NAN);
      out_.experiments.push_back(civ.experiments);
      out_.catastrophes.push_back(civ.catastrophes);
    }
    out_.unique_experiments = static_cast<long>(cache_.size());
    return std::move(out_);
  }

 private:
  Evaluation evaluate(const Genome& g) {
    std::string key = lab_.canonical(g);
    auto it = cache_.find(key);
    if (it != cache_.end()) return it->second;
    Evaluation e = lab_.evaluate(g);
    cache_.emplace(std::move(key), e);
    return e;
  }

  const Genome& pick_parent(Civ& civ) {
    if (!civ.catalog.empty() && rng_.chance(0.1)) {
      auto it = civ.catalog.begin();
      std::advance(it, rng_.below(static_cast<int>(civ.catalog.size())));
      return it->second.g;
    }
    // Torneio de 3: favorece bons resultados sem ignorar os demais.
    const Entry* winner = nullptr;
    for (int k = 0; k < 3; ++k) {
      const Entry& e = civ.archive[rng_.below(static_cast<int>(civ.archive.size()))];
      if (!winner || e.e.score > winner->e.score) winner = &e;
    }
    return winner->g;
  }

  void learn(int c, const Entry& entry, int year, bool by_contact) {
    Civ& civ = civs_[c];
    if (!entry.e.valid) return;

    if (!civ.has_best || entry.e.score > civ.best.e.score) {
      civ.best = entry;
      civ.has_best = true;
    }

    auto cat = civ.catalog.find(entry.e.label);
    if (cat == civ.catalog.end()) {
      civ.catalog.emplace(entry.e.label, entry);
      const bool first = known_labels_.insert(entry.e.label).second;
      out_.discoveries.push_back(
          {run_id_, year, c, entry.e.label, entry.e.score, first, by_contact, entry.g});
    } else if (entry.e.score > cat->second.e.score) {
      cat->second = entry;
    }

    for (const Entry& e : civ.archive)
      if (e.g == entry.g) return;
    if (static_cast<int>(civ.archive.size()) < cfg_.archive_size) {
      civ.archive.push_back(entry);
      return;
    }
    auto worst = std::min_element(civ.archive.begin(), civ.archive.end(),
                                  [](const Entry& a, const Entry& b) { return a.e.score < b.e.score; });
    if (entry.e.score > worst->e.score) *worst = entry;
  }

  void research_year(int c, int year) {
    Civ& civ = civs_[c];
    const Culture& cu = *civ.culture;
    double capacity =
        cu.speed * std::min(cfg_.max_experiments * cu.tools, cfg_.experiments0 * std::pow(1.0 + cfg_.growth, year));
    if (year < civ.dark_until) capacity *= 0.5;
    // Parte fracionária vira chance de um experimento a mais.
    int count = static_cast<int>(capacity);
    if (rng_.chance(capacity - count)) ++count;

    for (int i = 0; i < count; ++i) {
      Genome g;
      const double r = rng_.uniform();
      if (civ.archive.size() < 2 || r < civ.culture->explore) {
        g = lab_.random(rng_);
      } else if (r < civ.culture->explore + civ.culture->crossover) {
        g = lab_.crossover(pick_parent(civ), pick_parent(civ), rng_);
      } else {
        g = lab_.mutate(pick_parent(civ), rng_);
      }
      Entry entry{g, evaluate(g)};
      ++civ.experiments;
      learn(c, entry, year, false);
    }
  }

  void exchange(int c, int year) {
    Civ& civ = civs_[c];
    if (!civ.has_best || !rng_.chance(civ.culture->share)) return;
    std::vector<int> partners;
    for (int o = 0; o < static_cast<int>(civs_.size()); ++o)
      if (o != c && civs_[o].culture->share > 0) partners.push_back(o);
    if (partners.empty()) return;
    const int other = partners[rng_.below(static_cast<int>(partners.size()))];
    // Troca de mão dupla: o melhor resultado e uma classe aleatória do catálogo.
    for (auto [from, to] : {std::pair{c, other}, std::pair{other, c}}) {
      const Civ& src = civs_[from];
      if (!src.has_best) continue;
      learn(to, src.best, year, true);
      auto it = src.catalog.begin();
      std::advance(it, rng_.below(static_cast<int>(src.catalog.size())));
      learn(to, it->second, year, true);
    }
  }

  void maybe_catastrophe(int c, int year) {
    Civ& civ = civs_[c];
    if (!rng_.chance(cfg_.catastrophe_rate)) return;
    ++civ.catastrophes;
    civ.dark_until = year + cfg_.dark_age_years;
    std::vector<Entry> kept;
    for (auto& e : civ.archive)
      if (!rng_.chance(cfg_.catastrophe_loss)) kept.push_back(std::move(e));
    civ.archive = std::move(kept);
    // O catálogo (registro escrito) sobrevive parcialmente; o recorde só se estiver registrado.
    for (auto it = civ.catalog.begin(); it != civ.catalog.end();)
      it = rng_.chance(cfg_.catastrophe_loss * 0.5) ? civ.catalog.erase(it) : std::next(it);
    civ.has_best = false;
    for (const auto& e : civ.archive)
      if (!civ.has_best || e.e.score > civ.best.e.score) civ.best = e, civ.has_best = true;
    for (const auto& [label, e] : civ.catalog)
      if (!civ.has_best || e.e.score > civ.best.e.score) civ.best = e, civ.has_best = true;
  }

  const Lab& lab_;
  ResearchConfig cfg_;
  Rng rng_;
  int run_id_;
  std::vector<Civ> civs_;
  std::set<std::string> known_labels_;
  std::unordered_map<std::string, Evaluation> cache_;
  RunResult out_;
};

}  // namespace

RunResult run_research(const Lab& lab, const std::vector<Culture>& cultures, const ResearchConfig& cfg,
                       uint64_t seed, int run_id) {
  return Engine(lab, cultures, cfg, seed, run_id).run();
}
