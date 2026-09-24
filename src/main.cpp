// simu — civilizações simuladas pesquisando um universo regido pelo Modelo Padrão.
//
//   simu emergencia [--planetas N] [--saida DIR]
//   simu bioquimica [--planetas N] [--saida DIR]
//   simu pesquisa --lab brinquedo|quimica|materiais|nuclear [--anos Y] [--rodadas R]
//                 [--semente S] [--inicio E0] [--teto E] [--civs culturas|bioquimicas] [--saida DIR]
//   simu validar --lab ...
#include <omp.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <set>
#include <sstream>
#include <string>

#include "core/emergence.hpp"
#include "core/research.hpp"
#include "labs/labs.hpp"

namespace fs = std::filesystem;

namespace {

struct Args {
  std::string command;
  std::map<std::string, std::string> opts;
  std::string get(const std::string& k, const std::string& def) const {
    auto it = opts.find(k);
    return it == opts.end() ? def : it->second;
  }
};

Args parse(int argc, char** argv) {
  Args a;
  if (argc > 1) a.command = argv[1];
  for (int i = 2; i + 1 < argc; i += 2) {
    std::string k = argv[i];
    if (k.rfind("--", 0) == 0) k = k.substr(2);
    a.opts[k] = argv[i + 1];
  }
  return a;
}

std::unique_ptr<Lab> make_lab(const std::string& id, uint64_t seed) {
  if (id == "brinquedo") return make_toy_lab(seed);
  if (id == "quimica") return make_chemistry_lab();
  if (id == "materiais") return make_materials_lab();
  if (id == "nuclear") return make_nuclear_lab();
  throw std::runtime_error("laboratório desconhecido: " + id);
}

std::string csv_escape(std::string s) {
  std::replace(s.begin(), s.end(), '"', '\'');
  return "\"" + s + "\"";
}

double mean(const std::vector<double>& v) {
  return v.empty() ? NAN : std::accumulate(v.begin(), v.end(), 0.0) / v.size();
}
double stdev(const std::vector<double>& v) {
  if (v.size() < 2) return 0;
  const double m = mean(v);
  double s = 0;
  for (double x : v) s += (x - m) * (x - m);
  return std::sqrt(s / (v.size() - 1));
}
double quantile(std::vector<double> v, double q) {
  if (v.empty()) return NAN;
  std::sort(v.begin(), v.end());
  return v[static_cast<size_t>(q * (v.size() - 1))];
}

// ---------------------------------------------------------------- emergência
int cmd_emergence(const Args& args) {
  const long planets = std::stol(args.get("planetas", "20000000"));
  const fs::path dir = args.get("saida", "resultados/emergencia");
  fs::create_directories(dir);
  std::ofstream samples(dir / "amostras.csv");
  samples << "cenario,tempo_apos_formacao_gyr,tempo_cosmico_gyr\n";
  std::ofstream md(dir / "resumo.md");
  md << "# Camada A — Surgimento de civilizações científicas\n\n";
  md << "Monte Carlo com " << planets << " planetas rochosos na zona habitável por cenário. Formação dos sistemas "
     << "segue a taxa cósmica de formação estelar (Madau & Dickinson 2014) em ΛCDM (Planck 2018), ponderada pela "
     << "metalicidade; massa estelar 0,6–1,3 M☉ (IMF); janela habitável de 0,3 Gyr até 55% da vida da estrela.\n\n";
  md << "Referência — Terra: ciência começou ~4,54 Gyr após a formação do planeta, com o universo em 13,8 Gyr.\n\n";
  std::ostringstream table, steps;
  table << "| cenário | P(chegar à ciência) | mediana após formação (Gyr) | 10%–90% (Gyr) | "
           "idade do universo mediana (Gyr) | fração já possível até hoje |\n|---|---|---|---|---|---|\n";

  uint64_t seed = 1;
  for (const auto& sc : default_scenarios()) {
    const auto t0 = std::chrono::steady_clock::now();
    const EmergenceResult r = run_emergence(sc, planets, seed++);
    const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    const double p = static_cast<double>(r.successes) / r.planets;
    char buf[512];
    std::snprintf(buf, sizeof buf, "| %s | %.3g | %.2f | %.2f–%.2f | %.2f | %.2f |\n", sc.name.c_str(), p,
                  quantile(r.time_after_formation, 0.5), quantile(r.time_after_formation, 0.1),
                  quantile(r.time_after_formation, 0.9), quantile(r.cosmic_time, 0.5),
                  r.successes ? static_cast<double>(r.successes_before_today) / r.successes : 0.0);
    table << buf;
    std::cout << sc.name << ": P=" << p << " sucessos=" << r.successes << " (" << secs << " s)\n";
    for (size_t i = 0; i < r.time_after_formation.size() && i < 50000; ++i)
      samples << sc.name << "," << r.time_after_formation[i] << "," << r.cosmic_time[i] << "\n";

    // Cronologia mediana das transições nos planetas que chegaram lá.
    steps << "- **" << sc.name << "** (" << sc.description << "): ";
    for (size_t k = 0; k < sc.steps.size(); ++k)
      steps << sc.steps[k].name << " " << std::fixed << std::setprecision(2)
            << quantile(r.step_completion[k], 0.5) << " Gyr" << (k + 1 < sc.steps.size() ? " → " : "\n");
  }
  md << table.str() << "\n## Cronologia mediana (Gyr após a formação do planeta)\n\n" << steps.str();
  std::cout << "Resultados em " << dir << "\n";
  return 0;
}

// ---------------------------------------------------------------- bioquímicas
double median_or_nan(const std::vector<double>& v) { return v.empty() ? NAN : quantile(v, 0.5); }

// Civilizações de bioquímicas diferentes (cenário Terra, catálise adaptada): as 6 com mais chance de
// chegar à ciência. Velocidade = pensamento na temperatura mediana em que chegaram lá.
std::vector<Culture> biochemistry_civilizations() {
  const auto scenarios = default_scenarios();
  const EmergenceScenario& earth = scenarios[1];
  struct Cand { Culture c; double p; };
  std::vector<Cand> cands;
  uint64_t seed = 900;
  for (const auto& b : default_biochemistries()) {
    const auto r = run_bio_emergence(earth, b, 1.0, 400000, seed++);
    if (r.successes == 0) continue;
    const double t = median_or_nan(r.temperature);
    const BioFactors f = bio_factors(b, t, 1.0);
    // Cultura equilibrada, sem contato (estão em sistemas estelares diferentes).
    cands.push_back({{b.name, 0.15, 0.20, 0.0, f.thought, b.combustion ? 1.0 : 0.3},
                     static_cast<double>(r.successes) / r.planets});
  }
  std::sort(cands.begin(), cands.end(), [](const Cand& a, const Cand& b) { return a.p > b.p; });
  std::vector<Culture> out;
  for (size_t i = 0; i < cands.size() && i < 6; ++i) out.push_back(cands[i].c);
  return out;
}

int cmd_biochemistry(const Args& args) {
  const long planets = std::stol(args.get("planetas", "4000000"));
  const fs::path dir = args.get("saida", "resultados/bioquimica");
  fs::create_directories(dir);
  const auto biochems = default_biochemistries();
  const auto scenarios = default_scenarios();
  std::ofstream samples(dir / "amostras.csv");
  samples << "cenario,catalise,bioquimica,tempo_apos_formacao_gyr,tempo_cosmico_gyr,temperatura_k\n";
  std::ofstream md(dir / "resumo.md");
  md << "# Camada A — Surgimento com bioquímicas alternativas\n\n";
  md << "Monte Carlo com " << planets << " sistemas por bioquímica, cenário e hipótese de catálise. Cada planeta "
     << "tem órbita 0,05–50 UA, albedo 0,1–0,6 e efeito estufa ×1,00–1,35; a temperatura sobe com o brilho da "
     << "estrela, e a bioquímica só é viável com o solvente líquido. Os tempos de cada transição do cenário são "
     << "divididos pelos fatores químicos na temperatura do momento. **Química fixa**: Ea = 50 kJ/mol (frio = "
     << "lento). **Catálise adaptada**: Ea/RT igual ao da Terra (só o prefator kT/h muda).\n\n";

  md << "## Propriedades na temperatura típica do solvente\n\n";
  md << "| bioquímica | T típica (K) | reações, química fixa | reações, adaptada | estabilidade | versatilidade | "
        "polaridade | compatibilidade | fogo/metalurgia | nota |\n|---|---|---|---|---|---|---|---|---|---|\n";
  for (const auto& b : biochems) {
    const double t = b.solvent.name == "agua" ? 288.0 : 0.5 * (b.t_low() + b.t_high());
    const BioFactors f0 = bio_factors(b, t, 0.0), f1 = bio_factors(b, t, 1.0);
    char buf[512];
    std::snprintf(buf, sizeof buf, "| %s | %.0f | %.2g | %.2f | %.2f | %.2f | %.2f | %.2f | %s | %s |\n", b.name.c_str(),
                  t, f0.rate, f1.rate, f1.stability, f1.versatility, f1.polarity, b.compatibility,
                  b.combustion ? "sim" : "não", b.note.c_str());
    md << buf;
  }

  uint64_t seed = 77;
  for (const auto& sc : scenarios) {
    md << "\n## Cenário " << sc.name << "\n\n";
    md << "| bioquímica | planetas com ambiente viável | P(ciência), química fixa | P(ciência), adaptada | "
          "civilizações por civ. carbono/água (adaptada) | mediana até a ciência (Gyr, adaptada) | "
          "T mediana (K) | já possível até hoje |\n|---|---|---|---|---|---|---|---|\n";
    double reference = 0;
    std::vector<std::string> rows;
    for (const auto& b : biochems) {
      const auto r0 = run_bio_emergence(sc, b, 0.0, planets, seed++);
      const auto r1 = run_bio_emergence(sc, b, 1.0, planets, seed++);
      if (b.name == "Carbono/agua") reference = static_cast<double>(r1.successes);
      for (const auto* r : {&r0, &r1})
        for (size_t i = 0; i < r->time_after_formation.size() && i < 3000; ++i)
          samples << sc.name << "," << (r == &r0 ? "fixa" : "adaptada") << "," << b.name << ","
                  << r->time_after_formation[i] << "," << r->cosmic_time[i] << "," << r->temperature[i] << "\n";
      char buf[512];
      std::snprintf(buf, sizeof buf, "| %s | %.3f%% | %.2e | %.2e | %%s | %s | %s | %s |\n", b.name.c_str(),
                    100.0 * r1.habitable / r1.planets, static_cast<double>(r0.successes) / r0.planets,
                    static_cast<double>(r1.successes) / r1.planets,
                    r1.successes ? (std::to_string(median_or_nan(r1.time_after_formation)).substr(0, 5)).c_str() : "—",
                    r1.successes ? std::to_string(static_cast<int>(median_or_nan(r1.temperature))).c_str() : "—",
                    r1.successes ? (std::to_string(100 * r1.successes_before_today / r1.successes) + "%").c_str() : "—");
      rows.push_back(std::string(buf) + "\x1f" + std::to_string(r1.successes));
      std::cout << sc.name << " " << b.name << ": fixa " << r0.successes << ", adaptada " << r1.successes << "\n";
    }
    for (auto& row : rows) {
      const size_t cut = row.find('\x1f');
      const double n = std::stod(row.substr(cut + 1));
      char ratio[32];
      if (reference > 0) std::snprintf(ratio, sizeof ratio, "%.3g", n / reference);
      else std::snprintf(ratio, sizeof ratio, "—");
      char line[640];
      std::snprintf(line, sizeof line, row.substr(0, cut).c_str(), ratio);
      md << line;
    }
  }
  std::cout << "Resultados em " << dir << "\n";
  return 0;
}

// ---------------------------------------------------------------- pesquisa
std::string hex(const Genome& g) {
  static const char* digits = "0123456789abcdef";
  std::string s;
  for (uint8_t b : g) s += digits[b >> 4], s += digits[b & 15];
  return s;
}

std::string git_commit() {
  std::string out;
  if (FILE* p = popen("git rev-parse HEAD 2>/dev/null && git status --porcelain 2>/dev/null | head -1", "r")) {
    char buf[128];
    while (fgets(buf, sizeof buf, p)) out += buf;
    pclose(p);
  }
  const size_t nl = out.find('\n');
  std::string commit = out.substr(0, nl);
  const bool dirty = nl != std::string::npos && nl + 1 < out.size();
  return commit.empty() ? "desconhecido" : commit + (dirty ? " (com alterações locais)" : "");
}

// Catálogo reprodutível: parâmetros, semente, commit e, para cada classe, o genoma exato da primeira
// descoberta e do melhor exemplar, com a exportação completa do laboratório.
template <typename Items>
void write_catalog_json(const fs::path& path, const Lab& lab, const std::vector<Culture>& cultures,
                        const ResearchConfig& cfg, const Args& args, uint64_t seed, int runs, const Items& items) {
  std::ofstream js(path);
  js << "{\n  \"laboratorio\": \"" << lab.id() << "\",\n  \"commit\": \"" << git_commit() << "\",\n";
  js << "  \"comando\": \"simu pesquisa";
  for (const auto& [k, v] : args.opts) js << " --" << k << " " << v;
  js << "\",\n  \"semente\": " << seed << ", \"rodadas\": " << runs
     << ",\n  \"semente_de_cada_rodada\": \"semente * 1000003 + rodada\",\n";
  js << "  \"config\": {\"anos\": " << cfg.years << ", \"experimentos_ano0\": " << cfg.experiments0
     << ", \"crescimento\": " << cfg.growth << ", \"teto\": " << cfg.max_experiments
     << ", \"catastrofe_por_ano\": " << cfg.catastrophe_rate << ", \"perda_na_catastrofe\": " << cfg.catastrophe_loss
     << ", \"anos_de_trevas\": " << cfg.dark_age_years << ", \"memoria\": " << cfg.archive_size << "},\n";
  js << "  \"civilizacoes\": [";
  for (size_t c = 0; c < cultures.size(); ++c)
    js << (c ? ", " : "") << "{\"nome\": \"" << cultures[c].name << "\", \"explora\": " << cultures[c].explore
       << ", \"recombina\": " << cultures[c].crossover << ", \"contato\": " << cultures[c].share
       << ", \"velocidade\": " << cultures[c].speed << ", \"ferramentas\": " << cultures[c].tools << "}";
  js << "],\n  \"classes\": [\n";
  for (size_t i = 0; i < items.size(); ++i) {
    const auto& [label, it] = items[i];
    auto origin = [&](const auto& o) {
      std::ostringstream os;
      os << "{\"rodada\": " << o.run << ", \"civilizacao\": \"" << (o.civ >= 0 ? cultures[o.civ].name : "") << "\", \"ano\": "
         << o.year << ", \"como\": \"" << o.how << "\"}";
      return os.str();
    };
    js << "    {\"classe\": \"" << label << "\", \"melhor_score\": " << it.score << ", \"universos\": " << it.runs.size()
       << ",\n     \"melhor\": {\"origem\": " << origin(it.origin) << ", \"genoma_hex\": \"" << hex(it.genome)
       << "\", \"estrutura\": " << lab.export_json(it.genome) << "},\n     \"primeira\": {\"origem\": "
       << origin(it.first) << ", \"genoma_hex\": \"" << hex(it.first_genome) << "\", \"estrutura\": "
       << (it.first_genome.empty() ? std::string("{}") : lab.export_json(it.first_genome)) << "}}"
       << (i + 1 < items.size() ? ",\n" : "\n");
  }
  js << "  ]\n}\n";
}

int cmd_research(const Args& args) {
  const std::string lab_id = args.get("lab", "brinquedo");
  const uint64_t seed = std::stoull(args.get("semente", "42"));
  const int runs = std::stoi(args.get("rodadas", "8"));
  ResearchConfig cfg;
  cfg.years = std::stoi(args.get("anos", "300"));
  cfg.max_experiments = std::stod(args.get("teto", "1000"));
  cfg.experiments0 = std::stod(args.get("inicio", "10"));
  const auto lab = make_lab(lab_id, seed);
  const bool biochem_civs = args.get("civs", "culturas") == "bioquimicas";
  const auto cultures = biochem_civs ? biochemistry_civilizations() : default_cultures();
  const fs::path dir = args.get("saida", "resultados/" + lab_id + (biochem_civs ? "_bioquimicas" : ""));
  fs::create_directories(dir);

  std::cout << "Laboratório: " << lab->title() << "\n" << runs << " rodadas × " << cultures.size()
            << " civilizações × " << cfg.years << " anos, " << omp_get_max_threads() << " threads\n";
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<RunResult> results(runs);
#pragma omp parallel for schedule(dynamic)
  for (int r = 0; r < runs; ++r) results[r] = run_research(*lab, cultures, cfg, seed * 1000003 + r, r);
  const double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  // CSVs
  std::ofstream tl(dir / "linha_do_tempo.csv");
  tl << "rodada,ano,civilizacao,cultura,melhor,classes,experimentos\n";
  std::ofstream dc(dir / "descobertas.csv");
  dc << "rodada,ano,civilizacao,cultura,classe,score,primeira_do_universo,por_contato\n";
  long total_experiments = 0, unique = 0;
  double global_best = -1e300;
  for (const auto& res : results) {
    for (const auto& p : res.timeline)
      tl << p.run << "," << p.year << "," << p.civ << "," << cultures[p.civ].name << "," << p.best << ","
         << p.labels << "," << p.experiments << "\n";
    for (const auto& d : res.discoveries)
      dc << d.run << "," << d.year << "," << d.civ << "," << cultures[d.civ].name << "," << csv_escape(d.label) << ","
         << d.score << "," << d.global_first << "," << d.by_contact << "\n";
    for (long e : res.experiments) total_experiments += e;
    unique += res.unique_experiments;
    for (double s : res.best_score)
      if (!std::isnan(s)) global_best = std::max(global_best, s);
  }

  // Resumo por cultura
  std::ofstream md(dir / "resumo.md");
  md << "# " << lab->title() << "\n\n";
  md << "**Objetivo de pesquisa:** " << lab->objective() << ".\n\n";
  md << runs << " universos × " << cultures.size() << " civilizações × " << cfg.years << " anos de ciência; "
     << total_experiments << " experimentos (" << unique << " distintos) em " << std::fixed << std::setprecision(1)
     << secs << " s.\n\n";
  md << "## Validação contra o mundo real\n\n";
  for (const auto& line : lab->validate()) md << line << "\n";

  md << "\n## Desempenho por cultura de pesquisa\n\n";
  if (biochem_civs) {
    md << "Civilizações de bioquímicas diferentes (cenário Terra, catálise adaptada), cada uma no seu sistema "
          "estelar (sem contato). Velocidade = pensamento relativo ao humano; ferramentas = fração do teto "
          "instrumental (sem fogo/metalurgia: 30%).\n\n| civilização | velocidade | ferramentas |\n|---|---|---|\n";
    for (const auto& c : cultures) {
      char buf[160];
      std::snprintf(buf, sizeof buf, "| %s | %.2f | %.0f%% |\n", c.name.c_str(), c.speed, 100 * c.tools);
      md << buf;
    }
    md << "\n";
  }
  md << "| cultura | explora | recombina | contato | melhor score (média ± dp) | classes descobertas | "
        "descobertas inéditas no universo | ano médio em que atingiu o recorde global | catástrofes |\n";
  md << "|---|---|---|---|---|---|---|---|---|\n";
  for (size_t c = 0; c < cultures.size(); ++c) {
    std::vector<double> best, labels, reach, cats;
    int firsts = 0;
    for (const auto& res : results) {
      best.push_back(res.best_score[c]);
      cats.push_back(res.catastrophes[c]);
      std::set<std::string> distinct;
      for (const auto& d : res.discoveries)
        if (d.civ == static_cast<int>(c)) {
          distinct.insert(d.label);
          firsts += d.global_first;
        }
      labels.push_back(distinct.size());
      for (const auto& p : res.timeline)
        if (p.civ == static_cast<int>(c) && p.best >= global_best - 1e-9) {
          reach.push_back(p.year);
          break;
        }
    }
    char buf[512];
    std::snprintf(buf, sizeof buf, "| %s | %.2f | %.2f | %.2f | %.4g ± %.2g | %.1f | %d | %s (%zu/%d rodadas) | %.1f |\n",
                  cultures[c].name.c_str(), cultures[c].explore, cultures[c].crossover, cultures[c].share, mean(best),
                  stdev(best), mean(labels), firsts, reach.empty() ? "—" : std::to_string((int)mean(reach)).c_str(),
                  reach.size(), runs, mean(cats));
    md << buf;
  }

  // Catálogo de tudo o que foi descoberto
  struct Origin { int run = -1, civ = -1, year = -1; std::string how; };
  struct Item {
    int first_year = 1 << 30, runs_found = 0;
    double score = -1e300;
    Genome genome, first_genome;
    Origin origin, first;
    std::set<int> runs;
  };
  std::map<std::string, Item> catalog;
  for (const auto& res : results)
    for (const auto& d : res.discoveries) {
      Item& it = catalog[d.label];
      if (d.year < it.first_year) it.first_year = d.year, it.first_genome = d.genome, it.first = {d.run, d.civ, d.year, "descoberta"};
      it.runs.insert(d.run);
      if (d.score > it.score) it.score = d.score, it.genome = d.genome, it.origin = {d.run, d.civ, d.year, "descoberta"};
    }
  // Os melhores genomas finais também entram (podem ser melhores que o exemplo da 1ª descoberta).
  for (const auto& res : results)
    for (size_t c = 0; c < res.best_genome.size(); ++c) {
      if (res.best_genome[c].empty()) continue;
      const Evaluation e = lab->evaluate(res.best_genome[c]);
      Item& it = catalog[e.label];
      if (e.score > it.score)
        it.score = e.score, it.genome = res.best_genome[c],
        it.origin = {static_cast<int>(&res - results.data()), static_cast<int>(c), cfg.years, "melhor ao fim da rodada"};
    }
  std::vector<std::pair<std::string, Item>> items(catalog.begin(), catalog.end());
  std::sort(items.begin(), items.end(), [](auto& a, auto& b) { return a.second.score > b.second.score; });

  md << "\n## Catálogo de descobertas (" << items.size() << " classes; ordenadas pelo melhor score)\n\n";
  md << "| classe | melhor score | 1º ano (qualquer universo) | universos que acharam |\n|---|---|---|---|\n";
  for (const auto& [label, it] : items) {
    char buf[256];
    std::snprintf(buf, sizeof buf, "| %s | %.4g | %d | %zu/%d |\n", label.c_str(), it.score,
                  it.first_year == (1 << 30) ? -1 : it.first_year, it.runs.size(), runs);
    md << buf;
  }
  write_catalog_json(dir / "catalogo.json", *lab, cultures, cfg, args, seed, runs, items);

  md << "\n## Melhores exemplares\n\n";
  const int show = std::min<int>(8, items.size());
  for (int i = 0; i < show; ++i) {
    md << "**" << items[i].first << "** (score " << items[i].second.score << ")\n\n```\n"
       << lab->describe(items[i].second.genome) << "\n```\n\n";
  }

  md << "## História da ciência no universo 0 (primeiras descobertas)\n\n";
  for (const auto& d : results[0].discoveries)
    if (d.global_first)
      md << "- ano " << d.year << ": **" << cultures[d.civ].name << "** descobrem `" << d.label << "`\n";
  md.close();
  std::cout << "Concluído em " << secs << " s. Resultados em " << dir << "\n";
  return 0;
}

int cmd_validate(const Args& args) {
  const auto lab = make_lab(args.get("lab", "brinquedo"), 42);
  std::cout << lab->title() << "\n";
  for (const auto& line : lab->validate()) std::cout << line << "\n";
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  const Args args = parse(argc, argv);
  try {
    if (args.command == "emergencia") return cmd_emergence(args);
    if (args.command == "pesquisa") return cmd_research(args);
    if (args.command == "bioquimica") return cmd_biochemistry(args);
    if (args.command == "validar") return cmd_validate(args);
  } catch (const std::exception& e) {
    std::cerr << "erro: " << e.what() << "\n";
    return 1;
  }
  std::cerr << "uso: simu emergencia|pesquisa|validar [--lab brinquedo|quimica|materiais|nuclear] ...\n";
  return 1;
}
