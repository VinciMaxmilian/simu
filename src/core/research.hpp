#pragma once
#include <map>
#include <string>
#include <vector>

#include "core/lab.hpp"

// Cultura de pesquisa de uma civilização.
struct Culture {
  std::string name;
  double explore;    // fração de experimentos totalmente novos (curiosidade)
  double crossover;  // fração que combina duas ideias conhecidas
  double share;      // probabilidade anual de trocar conhecimento com outra civilização
  double speed = 1;  // velocidade de pensamento relativa à humana (bioquímica): multiplica os experimentos
  double tools = 1;  // fração do teto instrumental alcançável (sem fogo/metalurgia < 1)
};

std::vector<Culture> default_cultures();

struct ResearchConfig {
  int years = 300;                  // anos desde o início da ciência
  double experiments0 = 10;         // experimentos por ano no ano 0
  double growth = 0.02;             // crescimento anual da capacidade de pesquisa
  double max_experiments = 1000;    // teto de experimentos por ano
  double catastrophe_rate = 0.004;  // chance anual de uma idade das trevas
  double catastrophe_loss = 0.7;    // fração do acervo perdida numa catástrofe
  int dark_age_years = 25;          // anos com metade da capacidade após catástrofe
  int archive_size = 64;            // memória de trabalho (melhores resultados)
  int log_every = 5;
};

struct Discovery {
  int run, year, civ;
  std::string label;
  double score;
  bool global_first;  // ninguém havia descoberto antes
  bool by_contact;    // aprendido de outra civilização
  Genome genome;
};

struct TimelinePoint {
  int run, year, civ;
  double best;
  int labels;
  long experiments;
};

struct RunResult {
  std::vector<Discovery> discoveries;
  std::vector<TimelinePoint> timeline;
  std::vector<Genome> best_genome;  // por civilização
  std::vector<double> best_score;
  std::vector<long> experiments;
  std::vector<int> catastrophes;
  long unique_experiments = 0;
};

RunResult run_research(const Lab& lab, const std::vector<Culture>& cultures, const ResearchConfig& cfg,
                       uint64_t seed, int run_id);
