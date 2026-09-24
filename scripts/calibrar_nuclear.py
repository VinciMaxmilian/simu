"""Calibra o modelo de massas do laboratório nuclear (src/labs/nuclear.cpp).

B(Z,N) = a_v A − a_s A^(2/3) − a_c Z(Z−1)/A^(1/3) − a_a (N−Z)²/A + a_ss (N−Z)²/A^(4/3)
         + 12 δ/√A − S_ef        (emparelhamento fixo no valor usual de 12 MeV)
S      = C1 (F(Z)+F(N)) / (A/2)^(2/3) − C2 A^(1/3)          (Myers–Swiatecki)
S_ef   = S se S < 0, senão k·S  (núcleos no meio da camada se deformam e amortecem S > 0)

B é linear em (a_v, a_s, a_c, a_a, a_ss, C1, C2) para um padrão de sinais fixo de S:
resolvemos por mínimos quadrados iterando o padrão, para cada k numa grade.
Depois ajustamos log10 T_sf = s·B_f + i sobre os dados de fissão.
Uso: python3 scripts/calibrar_nuclear.py
"""
import csv, math, pathlib
import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
MAGIC_Z = [2, 8, 20, 28, 50, 82, 114, 164]
MAGIC_N = [2, 8, 20, 28, 50, 82, 126, 184, 258]


def shell_f(x, magic):
    prev = 0
    for m in magic:
        if x <= m:
            q = 0.6 * (m ** (5 / 3) - prev ** (5 / 3)) / (m - prev)
            return q * (x - prev) - 0.6 * (x ** (5 / 3) - prev ** (5 / 3))
        prev = m
    return 0.0


def shell_terms(z, n):
    a = z + n
    return (shell_f(z, MAGIC_Z) + shell_f(n, MAGIC_N)) / (a / 2) ** (2 / 3), a ** (1 / 3)


def features(z, n, p, k):
    """Linha de B(Z,N) nos 7 parâmetros, dado o vetor atual p (para decidir o sinal de S)."""
    a = z + n
    f, c = shell_terms(z, n)
    s = p[5] * f - p[6] * c
    w = 1.0 if s < 0 else k
    return np.array([a, -a ** (2 / 3), -z * (z - 1) / a ** (1 / 3), -(n - z) ** 2 / a, (n - z) ** 2 / a ** (4 / 3),
                     -w * f, w * c])


def pairing(z, n):
    delta = 1 if (z % 2 == 0 and n % 2 == 0) else -1 if (z % 2 == 1 and n % 2 == 1) else 0
    return 12.0 * delta / math.sqrt(z + n)


def binding(z, n, p, k):
    return features(z, n, p, k) @ p + pairing(z, n)


def load():
    rows = []
    with open(ROOT / "data/nuclear/calibracao.csv") as fh:
        for r in csv.DictReader(line for line in fh if not line.startswith("#")):
            rows.append((r["tipo"], int(r["Z"]), int(r["N"]), float(r["valor"])))
    return rows


def fit(rows, k):
    p = np.array([15.75, 17.8, 0.711, 23.7, 0.0, 5.8, 1.5])
    for _ in range(30):
        x, y, wts = [], [], []
        for t, z, n, v in rows:
            if t == "BA":
                x.append(features(z, n, p, k)); y.append(v * (z + n) - pairing(z, n)); wts.append(1.0)
            elif t == "QA":  # Qα = B(Z−2,N−2) + B(He4) − B(Z,N); peso maior: a meia-vida depende disso
                x.append(features(z - 2, n - 2, p, k) - features(z, n, p, k))
                y.append(v - 28.296 - pairing(z - 2, n - 2) + pairing(z, n)); wts.append(4.0)
        x, y, wts = np.array(x), np.array(y), np.array(wts)
        new = np.linalg.lstsq(x * wts[:, None], y * wts, rcond=None)[0]
        if np.allclose(new, p, atol=1e-9):
            break
        p = new
    res_b = [binding(z, n, p, k) / (z + n) - v for t, z, n, v in rows if t == "BA"]
    res_q = [binding(z - 2, n - 2, p, k) + 28.296 - binding(z, n, p, k) - v for t, z, n, v in rows if t == "QA"]
    return p, math.sqrt(np.mean(np.square(res_b))), math.sqrt(np.mean(np.square(res_q)))


def barrier(z, n, p, k):
    a = z + n
    es, ec = p[1] * a ** (2 / 3), p[2] * z * z / a ** (1 / 3)
    x = ec / (2 * es)
    bf = 0.38 * (0.75 - x) * es if x < 2 / 3 else 0.83 * (1 - x) ** 3 * es if x < 1 else 0.0
    f, c = shell_terms(z, n)
    s = p[5] * f - p[6] * c
    return bf - (s if s < 0 else k * s)


def main():
    rows = load()
    best = None
    for k in np.linspace(0, 1, 11):
        p, rb, rq = fit(rows, k)
        cost = rb * 10 + rq  # ~0.01 MeV em B/A ≈ importância de 0.1 MeV em Qα
        if best is None or cost < best[0]:
            best = (cost, k, p, rb, rq)
    _, k, p, rb, rq = best
    sf = [(barrier(z, n, p, k), v) for t, z, n, v in rows if t == "SF"]
    bfs, logs = np.array([b for b, _ in sf]), np.array([v for _, v in sf])
    slope, icpt = np.polyfit(bfs, logs, 1)
    rsf = math.sqrt(np.mean((slope * bfs + icpt - logs) ** 2))
    names = ["a_v", "a_s", "a_c", "a_a", "a_ss", "C1", "C2"]
    print(f"k_deformacao = {k:.2f}")
    for nm, v in zip(names, p):
        print(f"{nm:5s} = {v:.5f}")
    print(f"fissao: log10 T = {slope:.4f} * B_f + {icpt:.4f}")
    print(f"erro rms: B/A {rb * 1000:.1f} keV | Qalfa {rq:.3f} MeV | log10 T_sf {rsf:.2f}")


if __name__ == "__main__":
    main()
