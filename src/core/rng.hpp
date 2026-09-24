#pragma once
#include <cmath>
#include <cstdint>

// xoshiro256** — rápido, de boa qualidade e reprodutível a partir de uma semente.
class Rng {
 public:
  explicit Rng(uint64_t seed) {
    for (auto& word : s_) word = splitmix(seed);
  }

  uint64_t next() {
    const uint64_t result = rotl(s_[1] * 5, 7) * 9;
    const uint64_t t = s_[1] << 17;
    s_[2] ^= s_[0];
    s_[3] ^= s_[1];
    s_[1] ^= s_[2];
    s_[0] ^= s_[3];
    s_[2] ^= t;
    s_[3] = rotl(s_[3], 45);
    return result;
  }

  // Uniforme em [0, 1).
  double uniform() { return (next() >> 11) * 0x1.0p-53; }
  double uniform(double lo, double hi) { return lo + (hi - lo) * uniform(); }
  // Inteiro uniforme em [0, n).
  int below(int n) { return static_cast<int>(uniform() * n); }
  bool chance(double p) { return uniform() < p; }
  double exponential(double mean) { return -mean * std::log(1.0 - uniform()); }

 private:
  static uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
  static uint64_t splitmix(uint64_t& x) {
    uint64_t z = (x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
  }
  uint64_t s_[4];
};
