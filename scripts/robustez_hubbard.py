"""Robustez do resultado de Hubbard: sensibilidade a U, ao espaço ativo e à geometria real.

Para cada molécula (cálice de Clar = referência experimental; candidatos C49H21 e C51H25):
  • U/t ∈ {1.0, 1.2, 1.4, 1.6} no grafo plano;
  • espaço ativo do CAS: η+4 e η+6 orbitais;
  • saltos da geometria xTB real (Slater–Koster, parametrização de Moon–Koshino):
        t_ij = V_ppπ(d)·[n̂i·n̂j − (n̂i·b̂)(n̂j·b̂)] + V_ppσ(d)·(n̂i·b̂)(n̂j·b̂)
        V_ppπ(d) = −2.7 eV·exp(−(d − 1.42)/δ),  V_ppσ(d) = 0.48 eV·exp(−(d − 3.35)/δ),  δ = 0.453 Å
    em duas versões: só vizinhos ligados (o "t·cos φ") e todos os pares C–C até 4 Å (inclui o
    salto pelo espaço entre anéis sobrepostos de uma hélice).
O CAS é resolvido sem montar a matriz: H·C = Hα·C + C·Hβᵀ + U Σ_s Nα,s·C·Nβ,sᵀ, com Lanczos (eigsh).

Uso: python3 scripts/robustez_hubbard.py
Saída: resultados/candidatos/robustez_hubbard.json e .md
"""
import itertools, json, math, pathlib, sys, time

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import pipeline_quimica as pq  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "resultados" / "candidatos"
T0 = 2.7
DELTA = 0.184 * 2.46

MOLECULES = {
    "calice_de_clar": ("data/candidatos/C38H18_eta2.json", "calice_de_clar"),
    "C49H21": ("data/candidatos/C49H21_eta3_plano.json", "C49H21_eta3_plano"),
    "C51H25": ("data/candidatos/C51H25_eta3.json", "C51H25_eta3"),
}


# ------------------------------------------------------------------ saltos

def hopping_planar(struct):
    return -T0 * pq.adjacency(struct)


def xtb_carbons(folder, n_c):
    data = json.load(open(OUT / folder / "dados.json"))
    mult = min(data["xtb"], key=lambda k: data["xtb"][k]["energia_eV"])
    lines = [l for l in open(OUT / folder / f"geometria_xtb_mult{mult}.xyz").read().split("\n")[2:] if l.strip()]
    xyz = np.array([[float(v) for v in l.split()[1:4]] for l in lines])
    assert all(l.split()[0] == "C" for l in lines[:n_c])
    return xyz[:n_c]


def normals(xyz, adj):
    """Normal do plano local de cada carbono (vizinhos ligados), orientadas de forma coerente (BFS)."""
    n = len(xyz)
    nrm = np.zeros((n, 3))
    for i in range(n):
        pts = np.vstack([xyz[i], xyz[np.nonzero(adj[i])[0]]])
        _, _, vt = np.linalg.svd(pts - pts.mean(axis=0))
        nrm[i] = vt[2]
    seen, order = {0}, [0]
    for i in order:
        for j in np.nonzero(adj[i])[0]:
            if j not in seen:
                if nrm[i] @ nrm[j] < 0:
                    nrm[j] = -nrm[j]
                seen.add(j)
                order.append(j)
    return nrm


def hopping_geometry(struct, xyz, through_space):
    adj = pq.adjacency(struct)
    n = len(xyz)
    nrm = normals(xyz, adj)
    h = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = np.linalg.norm(xyz[j] - xyz[i])
            if not adj[i, j] and (not through_space or d > 4.0):
                continue
            b = (xyz[j] - xyz[i]) / d
            vpi = -T0 * math.exp(-(d - 1.42) / DELTA)
            vsig = 0.48 * math.exp(-(d - 3.35) / DELTA)
            ni_b, nj_b = nrm[i] @ b, nrm[j] @ b
            h[i, j] = h[j, i] = vpi * (nrm[i] @ nrm[j] - ni_b * nj_b) + vsig * ni_b * nj_b
    return h


# ------------------------------------------------------------------ CAS de Hubbard sem matriz

def cas_levels(h0, U, n_el, m, n_states=6):
    eps, phi = np.linalg.eigh(h0)
    lo = (n_el - m) // 2
    act = list(range(lo, lo + m))
    n_act = n_el - 2 * lo
    rho_core = (phi[:, :lo] ** 2).sum(axis=1)
    p = phi[:, act]
    h1 = np.diag(eps[act]) + U * (p.T * rho_core) @ p
    a_site = np.einsum("si,sj->sij", p, p)  # A_s[i,j] = φi(s)φj(s)

    def ops(k):
        strs = [sum(1 << i for i in c) for c in itertools.combinations(range(m), k)]
        index = {s: x for x, s in enumerate(strs)}
        e = np.zeros((m, m, len(strs), len(strs)))  # E_ij = a†_i a_j
        for x, s in enumerate(strs):
            for j in range(m):
                if not (s >> j) & 1:
                    continue
                sg = (-1) ** bin(s & ((1 << j) - 1)).count("1")
                s2 = s ^ (1 << j)
                for i in range(m):
                    if (s2 >> i) & 1:
                        continue
                    sg2 = sg * (-1) ** bin(s2 & ((1 << i) - 1)).count("1")
                    e[i, j, index[s2 | (1 << i)], x] += sg2
        return e

    levels = {}
    for n_a in range((n_act + 1) // 2, min(m, n_act) + 1):
        n_b = n_act - n_a
        if n_b < 0 or n_b > m:
            continue
        ea, eb = ops(n_a), ops(n_b)
        ha = np.einsum("ij,ijxy->xy", h1, ea)
        hb = np.einsum("ij,ijxy->xy", h1, eb)
        na = np.einsum("sij,ijxy->sxy", a_site, ea)
        nb = np.einsum("sij,ijxy->sxy", a_site, eb)
        da, db = ha.shape[0], hb.shape[0]

        def matvec(v):
            c = v.reshape(da, db)
            out = ha @ c + c @ hb.T
            out += U * np.einsum("sxy,yz,swz->xw", na, c, nb, optimize=True)
            return out.ravel()

        dim = da * db
        k = min(n_states, dim - 2) if dim > 3 else dim
        if dim <= 400:
            full = np.array([matvec(col) for col in np.eye(dim)]).T
            w = np.linalg.eigvalsh(full)[:n_states]
        else:
            w = np.sort(eigsh(LinearOperator((dim, dim), matvec=matvec, dtype=float), k=k, which="SA",
                              tol=1e-10, return_eigenvectors=False))
        levels[(n_a - n_b) / 2] = w
    # rótulo de S: um nível no setor Sz que não aparece no setor Sz+1 tem S = Sz
    out = []
    szs = sorted(levels)
    for idx, sz in enumerate(szs):
        higher = levels[szs[idx + 1]] if idx + 1 < len(szs) else []
        for e in levels[sz]:
            if not any(abs(e - e2) < 1e-6 for e2 in higher):
                out.append((float(e), sz))
    out.sort()
    e0 = out[0][0]
    return [(s, 1000 * (e - e0)) for e, s in out]


def first_excited_spin_gap(levels):
    s0 = levels[0][0]
    for s, e in levels[1:]:
        if s != s0:
            return s0, s, e
    return s0, None, None


def main():
    results = {}
    t_start = time.time()
    for name, (cand, folder) in MOLECULES.items():
        struct = json.load(open(ROOT / cand))["estrutura"]
        n_c = len(struct["carbonos"])
        eta = int((np.abs(np.linalg.eigvalsh(pq.adjacency(struct))) < 1e-8).sum())
        res = results[name] = {"n_C": n_c, "eta": eta}
        planar = hopping_planar(struct)
        # 1. varredura de U (grafo plano, espaço ativo η+4)
        res["U_scan"] = {}
        for u_t in (1.0, 1.2, 1.4, 1.6):
            lv = cas_levels(planar, u_t * T0, n_c, eta + 4)
            res["U_scan"][str(u_t)] = first_excited_spin_gap(lv)
            print(f"{name}: U/t={u_t} → {res['U_scan'][str(u_t)]}", flush=True)
        # 2. espaço ativo maior (U = 1.2 t)
        lv = cas_levels(planar, 1.2 * T0, n_c, eta + 6)
        res["CAS_maior"] = {"orbitais": eta + 6, "gap": first_excited_spin_gap(lv)}
        print(f"{name}: CAS({eta + 6}) → {res['CAS_maior']['gap']}", flush=True)
        # 3. geometria real (xTB)
        xyz = xtb_carbons(folder, n_c)
        for label, ts in (("geometria_vizinhos", False), ("geometria_com_espaco", True)):
            h = hopping_geometry(struct, xyz, ts)
            lv = cas_levels(h, 1.2 * T0, n_c, eta + 4)
            res[label] = first_excited_spin_gap(lv)
            t_bonds = [h[i, j] for i, j in struct["ligacoes"]]
            res[label + "_saltos_eV"] = [float(min(t_bonds)), float(max(t_bonds))]
            print(f"{name}: {label} → {res[label]}  (saltos ligados {min(t_bonds):.2f}…{max(t_bonds):.2f} eV)", flush=True)
    json.dump(results, open(OUT / "robustez_hubbard.json", "w"), indent=1)
    print(f"tempo total {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
