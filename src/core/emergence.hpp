#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "core/biochemistry.hpp"

// Camada A: quanto tempo leva, num universo com as leis do Modelo Padrão + ΛCDM,
// para um planeta chegar a uma civilização consciente que faz ciência?
//
// Modelo de "passos difíceis" (Carter 1983; Snyder-Beattie et al. 2021): a vida precisa
// completar uma sequência de transições, cada uma com tempo de espera exponencial,
// dentro da janela habitável do planeta.

struct EmergenceStep {
  std::string name;
  double mean_gyr;  // tempo médio de espera (bilhões de anos)
};

struct EmergenceScenario {
  std::string name;
  std::string description;
  std::vector<EmergenceStep> steps;
};

std::vector<EmergenceScenario> default_scenarios();

struct EmergenceResult {
  std::string scenario;
  long planets = 0;
  long successes = 0;
  long successes_before_today = 0;
  // Amostras (limitadas) dos sucessos, em Gyr.
  std::vector<double> time_after_formation;
  std::vector<double> cosmic_time;
  std::vector<std::vector<double>> step_completion;  // [passo][amostra]: fim do passo após a formação
};

EmergenceResult run_emergence(const EmergenceScenario& s, long planets, uint64_t seed);

// Surgimento para uma bioquímica específica. Cada planeta tem órbita (0.05–50 UA, log-uniforme),
// albedo e efeito estufa sorteados; a temperatura sobe com o brilho da estrela e a bioquímica só é
// viável enquanto o solvente está líquido. Os tempos de cada passo do cenário são divididos pelos
// fatores químicos (bio_factors) na temperatura do momento.
struct BioEmergenceResult {
  std::string biochemistry;
  long planets = 0, habitable = 0, successes = 0, successes_before_today = 0;
  std::vector<double> time_after_formation, cosmic_time, temperature;  // amostras dos sucessos
};

BioEmergenceResult run_bio_emergence(const EmergenceScenario& s, const Biochemistry& b, double catalysis,
                                     long planets, uint64_t seed);

// Idade do universo (Gyr) no redshift z, ΛCDM plano (Planck 2018).
double cosmic_time_gyr(double z);
