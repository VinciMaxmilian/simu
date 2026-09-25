#pragma once
#include <functional>
#include <utility>
#include <vector>

struct BenzenoidInfo {
  int rings = 0, carbons = 0, hydrogens = 0, n_a = 0, n_b = 0, eta = 0;
  double min_hh = 0;
  std::vector<std::pair<int, int>> cells;  // (q, r) normalizadas
};

struct EnumerationResult {
  std::vector<long> fixed, free;  // contagens por número de anéis (validação contra a OEIS)
};

// Visita cada benzenoide livre (a menos de rotação/reflexão) com até max_h anéis.
EnumerationResult enumerate_benzenoids(int max_h, const std::function<void(const BenzenoidInfo&)>& visit);
