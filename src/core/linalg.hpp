#pragma once
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <vector>

namespace linalg {

// Autovalores de uma matriz real simétrica n×n (row-major) pelo método de Jacobi.
// Retorna em ordem crescente. Adequado para matrizes pequenas (n ≲ 100).
inline std::vector<double> symmetric_eigenvalues(std::vector<double> a, int n) {
  auto at = [&](int i, int j) -> double& { return a[i * n + j]; };
  for (int sweep = 0; sweep < 100; ++sweep) {
    double off = 0.0;
    for (int i = 0; i < n; ++i)
      for (int j = i + 1; j < n; ++j) off += at(i, j) * at(i, j);
    if (off < 1e-22) break;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double apq = at(p, q);
        if (std::fabs(apq) < 1e-300) continue;
        const double theta = (at(q, q) - at(p, p)) / (2.0 * apq);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        for (int k = 0; k < n; ++k) {
          const double akp = at(k, p), akq = at(k, q);
          at(k, p) = c * akp - s * akq;
          at(k, q) = s * akp + c * akq;
        }
        for (int k = 0; k < n; ++k) {
          const double apk = at(p, k), aqk = at(q, k);
          at(p, k) = c * apk - s * aqk;
          at(q, k) = s * apk + c * aqk;
        }
      }
    }
  }
  std::vector<double> eig(n);
  for (int i = 0; i < n; ++i) eig[i] = at(i, i);
  std::sort(eig.begin(), eig.end());
  return eig;
}

// Autovalores de uma matriz hermitiana complexa n×n (row-major), Jacobi complexo.
// Cada rotação primeiro gira a fase de a_pq para torná-lo real e depois aplica a rotação real usual.
inline std::vector<double> hermitian_eigenvalues(std::vector<std::complex<double>> a, int n) {
  auto at = [&](int i, int j) -> std::complex<double>& { return a[i * n + j]; };
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = 0.0, diag = 0.0;
    for (int i = 0; i < n; ++i) {
      diag += std::norm(at(i, i));
      for (int j = i + 1; j < n; ++j) off += std::norm(at(i, j));
    }
    if (off <= 1e-30 * (diag + 1e-300)) break;
    for (int p = 0; p < n; ++p) {
      for (int q = p + 1; q < n; ++q) {
        const double mag = std::abs(at(p, q));
        if (mag < 1e-300) continue;
        const std::complex<double> w = std::conj(at(p, q)) / mag;  // e^{−iφ}
        const double theta = (at(q, q).real() - at(p, p).real()) / (2.0 * mag);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0);
        const double s = t * c;
        // A ← A·U com U_pp = c, U_qp = −s·w, U_pq = s, U_qq = c·w
        for (int k = 0; k < n; ++k) {
          const auto akp = at(k, p), akq = at(k, q);
          at(k, p) = c * akp - s * w * akq;
          at(k, q) = s * akp + c * w * akq;
        }
        // A ← Uᴴ·A
        for (int k = 0; k < n; ++k) {
          const auto apk = at(p, k), aqk = at(q, k);
          at(p, k) = c * apk - s * std::conj(w) * aqk;
          at(q, k) = s * apk + c * std::conj(w) * aqk;
        }
        at(p, q) = at(q, p) = 0.0;
      }
    }
  }
  std::vector<double> eig(n);
  for (int i = 0; i < n; ++i) eig[i] = at(i, i).real();
  std::sort(eig.begin(), eig.end());
  return eig;
}

// Caminho rápido para hermitianas pequenas (n ≤ 4), sem alocação: autovalores crescentes em out[0..n).
// n ≤ 2 em forma fechada; n = 3, 4 por Jacobi complexo em arrays fixos.
inline void small_hermitian_eigenvalues(const std::complex<double>* a, int n, double* out) {
  if (n == 1) {
    out[0] = a[0].real();
    return;
  }
  if (n == 2) {
    const double mean = 0.5 * (a[0].real() + a[3].real()), half = 0.5 * (a[0].real() - a[3].real());
    const double r = std::hypot(half, std::abs(a[1]));
    out[0] = mean - r, out[1] = mean + r;
    return;
  }
  std::complex<double> m[16];
  for (int i = 0; i < n * n; ++i) m[i] = a[i];
  auto at = [&](int i, int j) -> std::complex<double>& { return m[i * n + j]; };
  for (int sweep = 0; sweep < 50; ++sweep) {
    double off = 0.0, diag = 0.0;
    for (int i = 0; i < n; ++i) {
      diag += std::norm(at(i, i));
      for (int j = i + 1; j < n; ++j) off += std::norm(at(i, j));
    }
    if (off <= 1e-30 * (diag + 1e-300)) break;
    for (int p = 0; p < n; ++p)
      for (int q = p + 1; q < n; ++q) {
        const double mag = std::abs(at(p, q));
        if (mag < 1e-300) continue;
        const std::complex<double> w = std::conj(at(p, q)) / mag;
        const double theta = (at(q, q).real() - at(p, p).real()) / (2.0 * mag);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1.0));
        const double c = 1.0 / std::sqrt(t * t + 1.0), s = t * c;
        for (int k = 0; k < n; ++k) {
          const auto akp = at(k, p), akq = at(k, q);
          at(k, p) = c * akp - s * w * akq;
          at(k, q) = s * akp + c * w * akq;
        }
        for (int k = 0; k < n; ++k) {
          const auto apk = at(p, k), aqk = at(q, k);
          at(p, k) = c * apk - s * std::conj(w) * aqk;
          at(q, k) = s * apk + c * std::conj(w) * aqk;
        }
        at(p, q) = at(q, p) = 0.0;
      }
  }
  for (int i = 0; i < n; ++i) out[i] = at(i, i).real();
  std::sort(out, out + n);
}

// Determinante de uma matriz complexa n×n (eliminação gaussiana com pivotamento parcial).
inline std::complex<double> determinant(std::vector<std::complex<double>> a, int n) {
  std::complex<double> det = 1.0;
  for (int col = 0; col < n; ++col) {
    int pivot = col;
    for (int r = col + 1; r < n; ++r)
      if (std::abs(a[r * n + col]) > std::abs(a[pivot * n + col])) pivot = r;
    if (std::abs(a[pivot * n + col]) == 0.0) return 0.0;
    if (pivot != col) {
      for (int c = 0; c < n; ++c) std::swap(a[pivot * n + c], a[col * n + c]);
      det = -det;
    }
    det *= a[col * n + col];
    for (int r = col + 1; r < n; ++r) {
      const auto f = a[r * n + col] / a[col * n + col];
      for (int c = col; c < n; ++c) a[r * n + c] -= f * a[col * n + c];
    }
  }
  return det;
}

// Posto de uma matriz inteira (rows×cols, row-major) módulo um primo p.
inline int rank_mod_p(std::vector<int64_t> m, int rows, int cols, int64_t p) {
  auto mod_pow = [p](int64_t b, int64_t e) {
    int64_t r = 1;
    b %= p;
    while (e > 0) {
      if (e & 1) r = r * b % p;
      b = b * b % p;
      e >>= 1;
    }
    return r;
  };
  for (auto& x : m) x = ((x % p) + p) % p;
  int rank = 0;
  for (int col = 0; col < cols && rank < rows; ++col) {
    int pivot = -1;
    for (int r = rank; r < rows; ++r)
      if (m[r * cols + col] != 0) { pivot = r; break; }
    if (pivot < 0) continue;
    for (int c = 0; c < cols; ++c) std::swap(m[pivot * cols + c], m[rank * cols + c]);
    const int64_t inv = mod_pow(m[rank * cols + col], p - 2);
    for (int r = 0; r < rows; ++r) {
      if (r == rank || m[r * cols + col] == 0) continue;
      const int64_t f = m[r * cols + col] * inv % p;
      for (int c = col; c < cols; ++c)
        m[r * cols + c] = ((m[r * cols + c] - f * m[rank * cols + c]) % p + p) % p;
    }
    ++rank;
  }
  return rank;
}

// Posto sobre os racionais: o posto módulo p nunca excede o verdadeiro,
// então o máximo entre dois primos grandes é exato com probabilidade altíssima.
inline int rank_exact(const std::vector<int64_t>& m, int rows, int cols) {
  return std::max(rank_mod_p(m, rows, cols, 1000000007LL), rank_mod_p(m, rows, cols, 998244353LL));
}

}  // namespace linalg
