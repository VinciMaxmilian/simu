#pragma once
#include <string>
#include <vector>

// Bioquímicas possíveis: um elemento de esqueleto (as moléculas longas) + um solvente líquido.
//
// DADOS (físico-químicos, valores de referência):
//   energias médias de ligação (kJ/mol): C–C 346, C–O 358, C–H 411; Si–Si 222, Si–O 452, Si–H 318.
//   faixas líquidas a 1 atm (K) e constante dielétrica ε: água 273–373 (ε 80); amônia 195–240 (ε 22);
//   metano/etano 91–185 (ε ≈ 1.8); ácido sulfúrico 283–610 (ε ≈ 100); fluoreto de hidrogênio 190–293
//   (ε 84); nitrogênio líquido 63–77 (ε 1.43).
// HIPÓTESES (julgamento a partir da literatura; ajustáveis aqui):
//   abundância do solvente como oceano/lago, disponibilidade do elemento fora das rochas,
//   compatibilidade química esqueleto/solvente, e se há combustão/metalurgia (tecnologia).
//   Fontes: Bains 2004 (Astrobiology 4:137); Petkowski, Bains & Seager 2020 (Life 10:84, silício);
//   Sandström & Rahm 2020 (Sci. Adv., instabilidade dos azotossomos em Titã);
//   Seager et al. 2021 (PNAS, bases nucleicas estáveis em H2SO4 concentrado).

struct Backbone {
  std::string name;
  double e_chain;     // ligação X–X (kJ/mol): forma cadeias
  double e_oxide;     // ligação X–O (kJ/mol): concorrente (virar rocha/óxido)
  double e_weakest;   // ligação essencial mais fraca (kJ/mol): limita a estabilidade térmica
  double available;   // fração do elemento disponível fora de minerais (hipótese)
};

struct Solvent {
  std::string name;
  double t_melt, t_boil;  // K a 1 atm
  double dielectric;      // ε
  double abundance;       // probabilidade de um planeta ter oceano/lagos deste solvente (hipótese)
};

struct Biochemistry {
  std::string name;
  Backbone backbone;
  Solvent solvent;
  double compatibility;  // 0–1: o solvente preserva e permite reações do esqueleto (hipótese)
  bool combustion;       // há oxidante + combustível para fogo/metalurgia (hipótese)
  std::string note;

  // Faixa em que a superfície é habitável: do ponto de fusão até 80% do caminho à ebulição.
  double t_low() const { return solvent.t_melt; }
  double t_high() const { return solvent.t_melt + 0.8 * (solvent.t_boil - solvent.t_melt); }
};

std::vector<Biochemistry> default_biochemistries();

// Fatores relativos à Terra (carbono/água a 288 K = 1). catalysis ∈ [0,1]:
//   0 → química fixa (Ea = 50 kJ/mol, típica de enzimas);  1 → catálise adaptada (Ea/RT constante).
struct BioFactors {
  double rate;        // velocidade das reações (Arrhenius + prefator kT/h)
  double stability;   // sobrevivência das moléculas longas (ligação mais fraca vs kT)
  double versatility; // (E_cadeia/E_óxido)² relativo ao carbono
  double polarity;    // ε/(ε+10) relativo à água
  double thought;     // velocidade de pensamento = rate · stability
};

BioFactors bio_factors(const Biochemistry& b, double temperature, double catalysis);
