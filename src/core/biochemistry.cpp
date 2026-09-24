#include "core/biochemistry.hpp"

#include <algorithm>
#include <cmath>

namespace {

constexpr double kR = 8.314e-3;      // kJ/(mol·K)
constexpr double kEarthT = 288.0;    // K
constexpr double kEa = 50.0;         // kJ/mol, ativação típica de reações enzimáticas

const Backbone kCarbon{"carbono", 346, 358, 346, 1.0};
const Backbone kSilicon{"silicio", 222, 452, 222, 0.3};  // quase todo preso em silicatos

const Solvent kWater{"agua", 273.15, 373.15, 80.0, 1.0};
const Solvent kAmmonia{"amonia", 195.4, 239.8, 22.0, 0.3};
const Solvent kMethane{"metano/etano", 90.7, 184.6, 1.8, 0.5};
const Solvent kSulfuric{"acido sulfurico", 283.0, 610.0, 100.0, 0.05};
const Solvent kFluoride{"fluoreto de hidrogenio", 189.8, 292.7, 84.0, 0.01};
const Solvent kNitrogen{"nitrogenio liquido", 63.2, 77.4, 1.43, 0.1};

double polarity(double eps) { return eps / (eps + 10.0); }
double versatility(const Backbone& b) { return std::min(1.0, b.e_chain / b.e_oxide); }

}  // namespace

std::vector<Biochemistry> default_biochemistries() {
  return {
      {"Carbono/agua", kCarbon, kWater, 1.0, true, "a nossa; O2 de fotossíntese permite fogo e metalurgia"},
      {"Carbono/amonia", kCarbon, kAmmonia, 0.8, false, "química análoga à aquosa, mais fria; O2 reage com NH3"},
      {"Carbono/metano", kCarbon, kMethane, 0.3, false,
       "tipo Titã; solvente apolar dissolve pouco, membranas propostas (azotossomos) instáveis"},
      {"Carbono/acido sulfurico", kCarbon, kSulfuric, 0.3, false,
       "nuvens de Vênus; algumas moléculas orgânicas sobrevivem, sem superfície sólida"},
      {"Carbono/fluoreto", kCarbon, kFluoride, 0.3, false, "HF é bom solvente polar, mas flúor é raro"},
      {"Silicio/agua", kSilicon, kWater, 0.02, true, "silanos hidrolisam na água: vira sílica"},
      {"Silicio/metano", kSilicon, kMethane, 0.5, false, "silanos estáveis no frio e sem oxigênio"},
      {"Silicio/nitrogenio", kSilicon, kNitrogen, 0.5, false, "silanóis em N2 líquido (Bains 2004); muito frio"},
      {"Silicio/acido sulfurico", kSilicon, kSulfuric, 0.1, false, "quente o bastante para acelerar, mas corrosivo"},
  };
}

BioFactors bio_factors(const Biochemistry& b, double t, double catalysis) {
  BioFactors f;
  // Ea(T) interpola entre fixa (50 kJ/mol) e proporcional a T (mesmo Ea/RT da Terra).
  const double ea = (1 - catalysis) * kEa + catalysis * kEa * t / kEarthT;
  f.rate = (t / kEarthT) * std::exp(-ea / (kR * t) + kEa / (kR * kEarthT));
  // Estabilidade: moléculas longas precisam de ligações ≫ kT; limiar em E/RT ≈ 40.
  const double x = b.backbone.e_weakest / (kR * t);
  f.stability = 1.0 / (1.0 + std::exp(-(x - 40.0) / 4.0));
  f.versatility = std::pow(versatility(b.backbone) / versatility(kCarbon), 2.0);
  f.polarity = polarity(b.solvent.dielectric) / polarity(kWater.dielectric);
  f.thought = f.rate * f.stability;
  return f;
}
