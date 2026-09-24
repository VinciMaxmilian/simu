#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "core/rng.hpp"

// Um "experimento" é descrito por um genoma: bytes cujo significado é definido por cada laboratório.
using Genome = std::vector<uint8_t>;

struct Evaluation {
  bool valid = false;
  double score = -1e300;  // quanto maior, melhor (objetivo de pesquisa do laboratório)
  std::string label;      // classe qualitativa do resultado; a primeira ocorrência é uma "descoberta"
};

// Laboratório virtual: a física do universo simulado que as civilizações investigam.
// evaluate() precisa ser thread-safe (const, sem estado mutável).
class Lab {
 public:
  virtual ~Lab() = default;
  virtual std::string id() const = 0;
  virtual std::string title() const = 0;
  virtual std::string objective() const = 0;

  virtual Genome random(Rng& rng) const = 0;
  virtual Genome mutate(const Genome& g, Rng& rng) const = 0;
  virtual Genome crossover(const Genome& a, const Genome& b, Rng& rng) const = 0;
  virtual Evaluation evaluate(const Genome& g) const = 0;
  // Chave do cache: genomas com a mesma física devem ter a mesma chave (ex.: genes neutros,
  // translações). Padrão: os próprios bytes do genoma.
  virtual std::string canonical(const Genome& g) const { return std::string(g.begin(), g.end()); }
  virtual std::string describe(const Genome& g) const = 0;
  // Representação completa e legível por máquina (objeto JSON) para reconstruir o resultado depois.
  virtual std::string export_json(const Genome&) const { return "{}"; }

  // Confere o laboratório contra fatos conhecidos do mundo real (linhas em Markdown).
  virtual std::vector<std::string> validate() const = 0;
};
