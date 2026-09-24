#include "core/emergence.hpp"

#include <algorithm>
#include <cmath>
#include <omp.h>

#include "core/rng.hpp"

namespace {

constexpr double kH0 = 67.4;           // km/s/Mpc (Planck 2018)
constexpr double kOmegaM = 0.315;
constexpr double kHubbleTimeGyr = 977.8 / kH0;
constexpr double kToday = 13.8;        // Gyr
constexpr int kMaxSamples = 200000;

// Taxa cósmica de formação estelar (Madau & Dickinson 2014), em unidades arbitrárias.
double star_formation_rate(double z) {
  return std::pow(1.0 + z, 2.7) / (1.0 + std::pow((1.0 + z) / 2.9, 5.6));
}

// Fração de estrelas com metais suficientes para planetas rochosos (aproximação suave).
double metal_fraction(double t_gyr) { return 1.0 - std::exp(-t_gyr / 2.0); }

// Tabela: tempo cósmico → peso de formação de sistemas planetários habitáveis.
struct FormationSampler {
  std::vector<double> t, cdf;
  FormationSampler() {
    const int n = 2000;
    for (int i = 0; i <= n; ++i) {
      const double z = 20.0 * std::pow(static_cast<double>(n - i) / n, 2.0);  // mais pontos em z baixo
      t.push_back(cosmic_time_gyr(z));
    }
    cdf.assign(t.size(), 0.0);
    for (size_t i = 1; i < t.size(); ++i) {
      auto w = [&](size_t k) {
        const double z = 20.0 * std::pow(static_cast<double>(t.size() - 1 - k) / (t.size() - 1), 2.0);
        return star_formation_rate(z) * metal_fraction(t[k]);
      };
      cdf[i] = cdf[i - 1] + 0.5 * (w(i) + w(i - 1)) * (t[i] - t[i - 1]);
    }
    for (auto& c : cdf) c /= cdf.back();
  }
  double sample(Rng& rng) const {
    const double u = rng.uniform();
    const size_t i = std::lower_bound(cdf.begin(), cdf.end(), u) - cdf.begin();
    if (i == 0) return t[0];
    const double f = (u - cdf[i - 1]) / std::max(1e-300, cdf[i] - cdf[i - 1]);
    return t[i - 1] + f * (t[i] - t[i - 1]);
  }
};

// Massa estelar (massas solares) com IMF de Salpeter/Kroupa (inclinação −2.3), 0.6–1.3 M☉.
double sample_star_mass(Rng& rng) {
  const double a = 1.0 - 2.3, lo = std::pow(0.6, a), hi = std::pow(1.3, a);
  return std::pow(lo + (hi - lo) * rng.uniform(), 1.0 / a);
}

}  // namespace

double cosmic_time_gyr(double z) {
  // t(z) = ∫_z^∞ dz' / ((1+z') H(z')), resolvido analiticamente para ΛCDM plano.
  const double omega_l = 1.0 - kOmegaM;
  const double x = std::sqrt(omega_l / kOmegaM) * std::pow(1.0 + z, -1.5);
  return kHubbleTimeGyr * 2.0 / (3.0 * std::sqrt(omega_l)) * std::asinh(x);
}

std::vector<EmergenceScenario> default_scenarios() {
  const double sci = 0.0003;  // da inteligência à ciência: ~300 mil anos (Homo sapiens → método científico)
  return {
      {"Rapido", "vida complexa e inteligência são fáceis; cada transição leva centenas de milhões de anos",
       {{"abiogenese", 0.2}, {"eucariontes", 0.8}, {"multicelularidade complexa", 0.6},
        {"inteligencia/linguagem", 0.4}, {"ciencia", sci}}},
      {"Terra", "tempos médios iguais aos intervalos observados na Terra",
       {{"abiogenese", 0.5}, {"eucariontes", 1.8}, {"multicelularidade complexa", 1.2},
        {"inteligencia/linguagem", 0.55}, {"ciencia", sci}}},
      {"PassosDificeis", "transições raras (tempo médio ≫ janela habitável), como sugere o modelo de passos difíceis",
       {{"abiogenese", 10.0}, {"eucariontes", 10.0}, {"multicelularidade complexa", 10.0},
        {"inteligencia/linguagem", 10.0}, {"ciencia", sci}}},
  };
}

EmergenceResult run_emergence(const EmergenceScenario& s, long planets, uint64_t seed) {
  static const FormationSampler sampler;
  EmergenceResult res;
  res.scenario = s.name;
  res.planets = planets;
  res.step_completion.resize(s.steps.size());
  const int threads = omp_get_max_threads();

#pragma omp parallel num_threads(threads)
  {
    const int tid = omp_get_thread_num();
    Rng rng(seed * 7919 + tid);
    EmergenceResult local;
    local.step_completion.resize(s.steps.size());
    std::vector<double> ends(s.steps.size());
#pragma omp for schedule(static)
    for (long p = 0; p < planets; ++p) {
      const double born = sampler.sample(rng);
      const double mass = sample_star_mass(rng);
      const double lifetime = 10.0 * std::pow(mass, -2.5);
      // Janela habitável: do resfriamento da crosta (~0.3 Gyr) até o brilho crescente da
      // estrela esterilizar o planeta (~55% da vida na sequência principal; Terra: ~5.5 Gyr).
      const double window_end = 0.55 * lifetime;
      double t = 0.3;
      bool ok = true;
      for (size_t k = 0; k < s.steps.size(); ++k) {
        t += rng.exponential(s.steps[k].mean_gyr);
        ends[k] = t;
        if (t > window_end) { ok = false; break; }
      }
      if (!ok) continue;
      ++local.successes;
      if (born + t <= kToday) ++local.successes_before_today;
      if (static_cast<long>(local.time_after_formation.size()) < kMaxSamples / threads) {
        local.time_after_formation.push_back(t);
        local.cosmic_time.push_back(born + t);
        for (size_t k = 0; k < ends.size(); ++k) local.step_completion[k].push_back(ends[k]);
      }
    }
#pragma omp critical
    {
      res.successes += local.successes;
      res.successes_before_today += local.successes_before_today;
      auto append = [](std::vector<double>& a, const std::vector<double>& b) { a.insert(a.end(), b.begin(), b.end()); };
      append(res.time_after_formation, local.time_after_formation);
      append(res.cosmic_time, local.cosmic_time);
      for (size_t k = 0; k < ends.size(); ++k) append(res.step_completion[k], local.step_completion[k]);
    }
  }
  return res;
}
