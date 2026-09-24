"""Pipeline de verificação de candidatos da química: Hückel → Hubbard → xTB → DFT.

Parte do grafo EXATO exportado pelo simulador (data/candidatos/*.json ou um catálogo) e responde:
  1. Hubbard: o magnetismo sobrevive à repulsão elétron-elétron? (campo médio + CAS exato)
  2. xTB (GFN2): geometria otimizada em cada estado de spin, gap e frequências (mínimo estável?)
  3. DFT spin-polarizada (UB3LYP/6-31G): energias dos estados de spin, <S²> e densidade de spin

Uso:
  python3 scripts/pipeline_quimica.py data/candidatos/C51H25_eta3.json [--sem-dft] [--sem-freq]
  python3 scripts/pipeline_quimica.py --referencias      # valida em moléculas com resultado conhecido
Saída: resultados/candidatos/<nome>/ (relatorio.md, dados.json, geometrias .xyz, densidade de spin .png)
"""
import argparse, itertools, json, math, pathlib, sys, time

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "resultados" / "candidatos"
T_HOP = 2.7         # eV, salto entre carbonos vizinhos (grafeno)
U_HUB = 1.2 * T_HOP  # eV, repulsão local típica em nanografenos (Mishra et al. 2020; Ortiz et al.)
CC, CH = 1.42, 1.09  # Å

# ------------------------------------------------------------------ estrutura

# Mesma geometria do simulador (src/labs/chemistry.cpp): hexágonos "pointy-top" em coordenadas axiais.
CORNERS = [(0, 2), (1, 1), (1, -1), (0, -2), (-1, -1), (-1, 1)]


def structure_from_cells(cells):
    atoms, pos, bonds = {}, [], set()
    for q, r in cells:
        cx, cy = 2 * q + r, 3 * r
        ids = []
        for dx, dy in CORNERS:
            key = (cx + dx, cy + dy)
            if key not in atoms:
                atoms[key] = len(pos)
                pos.append(key)
            ids.append(atoms[key])
        for k in range(6):
            bonds.add(tuple(sorted((ids[k], ids[(k + 1) % 6]))))
    carbons = [{"id": i, "x_A": x * math.sqrt(3) / 2 * CC, "y_A": y * 0.5 * CC,
                "subrede": "A" if (y % 3) == 2 else "B"} for i, (x, y) in enumerate(pos)]
    return {"anel_coordenadas_axiais_q_r": [list(c) for c in cells], "carbonos": carbons,
            "ligacoes": [list(b) for b in sorted(bonds)]}


def adjacency(struct):
    n = len(struct["carbonos"])
    a = np.zeros((n, n))
    for i, j in struct["ligacoes"]:
        a[i, j] = a[j, i] = 1
    return a


def geometry(struct):
    """Carbonos no plano + H nas bordas (carbonos com 2 vizinhos), apontando para fora."""
    c = np.array([[a["x_A"], a["y_A"], 0.0] for a in struct["carbonos"]])
    a = adjacency(struct)
    symbols, coords = ["C"] * len(c), list(c)
    for i in range(len(c)):
        nb = np.nonzero(a[i])[0]
        if len(nb) == 2:
            v = c[i] - c[nb].mean(axis=0)
            coords.append(c[i] + CH * v / np.linalg.norm(v))
            symbols.append("H")
    return symbols, np.array(coords)


def write_xyz(path, symbols, coords, comment=""):
    with open(path, "w") as fh:
        fh.write(f"{len(symbols)}\n{comment}\n")
        for s, (x, y, z) in zip(symbols, coords):
            fh.write(f"{s} {x:.6f} {y:.6f} {z:.6f}\n")


# ------------------------------------------------------------------ Hubbard

def hubbard_mf(a, n_up, n_dn, guess, iters=3000, mix=0.3, tol=1e-9):
    """Hubbard em campo médio (UHF). Retorna energia, densidades por spin, convergiu?"""
    h0 = -T_HOP * a
    up, dn = guess
    for _ in range(iters):
        eu, vu = np.linalg.eigh(h0 + U_HUB * np.diag(dn))
        ed, vd = np.linalg.eigh(h0 + U_HUB * np.diag(up))
        new_up = (vu[:, :n_up] ** 2).sum(axis=1)
        new_dn = (vd[:, :n_dn] ** 2).sum(axis=1)
        delta = max(np.abs(new_up - up).max(), np.abs(new_dn - dn).max())
        up, dn = (1 - mix) * up + mix * new_up, (1 - mix) * dn + mix * new_dn
        if delta < tol:
            break
    energy = eu[:n_up].sum() + ed[:n_dn].sum() - U_HUB * (up * dn).sum()
    return energy, up, dn, delta < 1e-6


def hubbard_mf_scan(struct, n_electrons, seed=0):
    """Menor energia de campo médio para cada Sz, a partir de várias condições iniciais."""
    a = adjacency(struct)
    n = len(a)
    sub = np.array([1.0 if c["subrede"] == "A" else -1.0 for c in struct["carbonos"]])
    rng = np.random.default_rng(seed)
    results = {}
    max_sz2 = min(8, n_electrons)  # 2·Sz até 4 spins alinhados (suficiente aqui)
    for two_sz in range(n_electrons % 2, max_sz2 + 1, 2):
        n_up, n_dn = (n_electrons + two_sz) // 2, (n_electrons - two_sz) // 2
        base = np.full(n, n_electrons / (2 * n))
        guesses = [(base + 0.3 * sub, base - 0.3 * sub), (base - 0.3 * sub, base + 0.3 * sub),
                   (base + 0.1, base - 0.1)]
        guesses += [(base + 0.3 * r, base - 0.3 * r) for r in rng.normal(size=(5, n))]
        best = None
        for g in guesses:
            e, up, dn, ok = hubbard_mf(a, n_up, n_dn, (np.clip(g[0], 0, 1), np.clip(g[1], 0, 1)))
            if best is None or e < best[0] - 1e-9:
                best = (e, up, dn, ok)
        results[two_sz / 2] = {"energia_eV": best[0], "densidade_spin": (best[1] - best[2]).tolist(),
                               "momento_total": float((best[1] - best[2]).sum()),
                               "momento_local_max": float(np.abs(best[1] - best[2]).max()), "convergiu": bool(best[3])}
    return results


def hubbard_cas(struct, n_electrons, n_active):
    """Diagonalização exata do Hubbard num espaço ativo de orbitais de Hückel em torno de E = 0.
    Caroço duplamente ocupado entra como campo médio (U·densidade do caroço). Retorna energias por S."""
    a = adjacency(struct)
    eps, phi = np.linalg.eigh(-T_HOP * a)
    order = np.argsort(np.abs(eps), kind="stable")[:n_active]
    act = np.sort(order)
    core = [i for i in range(len(eps)) if i < act.min()]
    assert all(i not in core for i in act)
    n_act_el = n_electrons - 2 * len(core)
    rho_core = (phi[:, core] ** 2).sum(axis=1)
    p = phi[:, act]
    h1 = np.diag(eps[act]) + U_HUB * (p.T * rho_core) @ p          # h_ij = ε_i δ_ij + U Σ_s φiφj ρ_core
    v = U_HUB * np.einsum("si,sj,sk,sl->ijkl", p, p, p, p)          # (ij|kl) do termo n↑n↓
    m = n_active

    def strings(k):
        return [sum(1 << i for i in c) for c in itertools.combinations(range(m), k)]

    def apply(i, j, s):  # a†_i a_j |s>: (sinal, novo) ou None
        if not (s >> j) & 1:
            return None
        sign = (-1) ** bin(s & ((1 << j) - 1)).count("1")
        s2 = s ^ (1 << j)
        if (s2 >> i) & 1:
            return None
        sign *= (-1) ** bin(s2 & ((1 << i) - 1)).count("1")
        return sign, s2 | (1 << i)

    energies = {}
    for n_a in range((n_act_el + 1) // 2, n_act_el + 1):
        n_b = n_act_el - n_a
        if n_a > m or n_b < 0:
            continue
        sa, sb = strings(n_a), strings(n_b)
        ia, ib = {s: k for k, s in enumerate(sa)}, {s: k for k, s in enumerate(sb)}
        dim = len(sa) * len(sb)
        # Excitações de um corpo (i←j) para cada string, pré-calculadas.
        ex_a = [[(i, j, *r) for i in range(m) for j in range(m) if (r := apply(i, j, s))] for s in sa]
        ex_b = [[(i, j, *r) for i in range(m) for j in range(m) if (r := apply(i, j, s))] for s in sb]
        h = np.zeros((dim, dim))
        for x, s_a in enumerate(sa):
            for y, s_b in enumerate(sb):
                col = x * len(sb) + y
                for i, j, sg, t in ex_a[x]:
                    h[ia[t] * len(sb) + y, col] += h1[i, j] * sg
                for i, j, sg, t in ex_b[y]:
                    h[x * len(sb) + ib[t], col] += h1[i, j] * sg
                for i, j, sg1, t1 in ex_a[x]:
                    for k, l, sg2, t2 in ex_b[y]:
                        h[ia[t1] * len(sb) + ib[t2], col] += v[i, j, k, l] * sg1 * sg2
        w, vec = np.linalg.eigh(h)
        two_sz = n_a - n_b
        # S de cada autoestado: o menor Sz em que aparece define S; guardamos estados com S = Sz aqui
        for e in w:
            energies.setdefault(two_sz / 2, []).append(float(e))
    # Um multipleto de spin S aparece em todo setor com Sz ≤ S. Energia do menor S possível:
    levels = []
    sz_values = sorted(energies)
    for idx, sz in enumerate(sz_values):
        higher = energies.get(sz_values[idx + 1], []) if idx + 1 < len(sz_values) else []
        for e in energies[sz]:
            # estados deste setor que não existem no setor Sz+1 têm S = Sz
            if not any(abs(e - e2) < 1e-6 for e2 in higher):
                levels.append((e, sz))
    levels.sort()
    e0 = levels[0][0]
    return {"orbitais_ativos": int(n_active), "eletrons_ativos": int(n_act_el),
            "niveis": [{"S": s, "energia_rel_meV": 1000 * (e - e0)} for e, s in levels[:8]]}


# ------------------------------------------------------------------ xTB e DFT

def xtb_optimize(symbols, coords, multiplicity, fmax=0.03, freq=False):
    from ase import Atoms
    from ase.optimize import BFGS
    from tblite.ase import TBLite
    atoms = Atoms(symbols, positions=coords)
    atoms.calc = TBLite(method="GFN2-xTB", multiplicity=multiplicity, verbosity=0)
    # Pequena perturbação fora do plano: deixa a molécula dobrar se o plano não for mínimo.
    atoms.positions[:, 2] += np.random.default_rng(1).normal(scale=0.02, size=len(atoms))
    BFGS(atoms, logfile=None).run(fmax=fmax, steps=2000)
    res = {"multiplicidade": multiplicity, "energia_eV": float(atoms.get_potential_energy()),
           "desvio_do_plano_A": float(np.abs(atoms.positions[:, 2] - atoms.positions[:, 2].mean()).max())}
    cc = [np.linalg.norm(atoms.positions[i] - atoms.positions[j]) for i in range(len(atoms)) for j in range(i)
          if symbols[i] == symbols[j] == "C" and np.linalg.norm(atoms.positions[i] - atoms.positions[j]) < 1.7]
    res["ligacao_CC_A"] = [float(min(cc)), float(max(cc))]
    calc = atoms.calc
    try:
        eps = calc.results.get("orbital_energies")
        occ = calc.results.get("orbital_occupations")
    except Exception:
        eps = occ = None
    if freq:
        from ase.vibrations import Vibrations
        import shutil, tempfile
        tmp = tempfile.mkdtemp()
        vib = Vibrations(atoms, name=f"{tmp}/vib", delta=0.005)
        vib.run()
        f = vib.get_frequencies()
        shutil.rmtree(tmp)
        real = np.real(f)[np.abs(np.imag(f)) < 1e-6]
        imag = np.imag(f)[np.abs(np.imag(f)) > 1e-6]
        res["frequencias"] = {"menor_real_cm1": float(np.sort(real[real > 1])[0]) if len(real[real > 1]) else None,
                              "imaginarias_cm1": [float(x) for x in imag if abs(x) > 30],
                              "n_modos": int(len(f))}
    return res, atoms.positions.copy()


def dft_single_point(symbols, coords, spin2s, basis="6-31g", xc="b3lyp", breaksym=False):
    from pyscf import dft, gto, lib
    lib.num_threads(4)
    mol = gto.M(atom=[(s, tuple(c)) for s, c in zip(symbols, coords)], basis=basis, spin=spin2s, verbose=0)
    mf = dft.UKS(mol).density_fit()
    mf.xc = xc
    mf.max_cycle = 200
    mf.conv_tol = 1e-7
    if breaksym:
        mf.init_guess_breaksym = True
    mf.kernel()
    s2, mult = mf.spin_square()
    pop, _ = mf.mulliken_pop(verbose=0)
    dm_a, dm_b = mf.make_rdm1()
    spin_pop = mf.mulliken_pop(dm=dm_a - dm_b, verbose=0)[1] if False else None
    # densidade de spin de Mulliken por átomo
    s = mol.intor_symmetric("int1e_ovlp")
    ao_spin = np.einsum("ij,ji->i", dm_a - dm_b, s)
    per_atom = np.zeros(mol.natm)
    for i, (ia, *_rest) in enumerate(mol.ao_labels(fmt=False)):
        per_atom[ia] += ao_spin[i]
    mo_e = np.concatenate(mf.mo_energy)
    occ = np.concatenate(mf.mo_occ)
    homo, lumo = mo_e[occ > 0.5].max(), mo_e[occ < 0.5].min()
    return {"2S": spin2s, "quebra_de_simetria": breaksym, "energia_Ha": float(mf.e_tot), "convergiu": bool(mf.converged),
            "S2": float(s2), "gap_somo_lumo_eV": float(27.2114 * (lumo - homo)),
            "densidade_spin_atomos": per_atom.tolist()}


# ------------------------------------------------------------------ relatório

def spin_density_plot(path, struct, spin_c, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    xy = np.array([[c["x_A"], c["y_A"]] for c in struct["carbonos"]])
    fig, ax = plt.subplots(figsize=(6, 5))
    for i, j in struct["ligacoes"]:
        ax.plot(*xy[[i, j]].T, color="#9aa1b0", lw=1, zorder=1)
    m = max(1e-9, np.abs(spin_c).max())
    ax.scatter(xy[:, 0], xy[:, 1], s=900 * np.abs(spin_c) / m + 4,
               c=np.where(np.array(spin_c) >= 0, "#2a78d6", "#eb6834"), alpha=0.8, zorder=2, edgecolors="none")
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_title(title, fontsize=10, loc="left")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def run(name, struct, args, expected=None):
    out = OUT / name
    out.mkdir(parents=True, exist_ok=True)
    a = adjacency(struct)
    n_c = len(a)
    subs = [c["subrede"] for c in struct["carbonos"]]
    n_a, n_b = subs.count("A"), subs.count("B")
    eig = np.linalg.eigvalsh(a)
    eta = int((np.abs(eig) < 1e-8).sum())
    symbols, coords = geometry(struct)
    data = {"nome": name, "formula": f"C{n_c}H{symbols.count('H')}", "n_A": n_a, "n_B": n_b, "eta_huckel": eta,
            "S_lieb": abs(n_a - n_b) / 2, "parametros": {"t_eV": T_HOP, "U_eV": U_HUB}}
    if expected:
        data["esperado_literatura"] = expected
    t0 = time.time()
    print(f"[{name}] {data['formula']}  η={eta}  N_A/N_B={n_a}/{n_b}  Lieb: S={data['S_lieb']}", flush=True)

    # 1. Hubbard
    mf = hubbard_mf_scan(struct, n_c)
    e_min = min(v["energia_eV"] for v in mf.values())
    data["hubbard_campo_medio"] = {str(k): {kk: vv for kk, vv in v.items() if kk != "densidade_spin"}
                                   | {"energia_rel_meV": 1000 * (v["energia_eV"] - e_min)} for k, v in mf.items()}
    gs_sz = min(mf, key=lambda k: mf[k]["energia_eV"])
    data["hubbard_campo_medio_Sz_fundamental"] = gs_sz
    n_active = min(8, eta + 4) if eta else 4
    data["hubbard_cas"] = hubbard_cas(struct, n_c, n_active)
    print(f"  Hubbard MF: Sz fundamental={gs_sz}, momento local máx={mf[gs_sz]['momento_local_max']:.3f}", flush=True)
    print(f"  Hubbard CAS({data['hubbard_cas']['eletrons_ativos']},{n_active}): "
          + ", ".join(f"S={l['S']}:{l['energia_rel_meV']:.1f} meV" for l in data["hubbard_cas"]["niveis"][:4]), flush=True)
    spin_density_plot(out / "densidade_spin_hubbard.png", struct, mf[gs_sz]["densidade_spin"],
                      f"{data['formula']} — densidade de spin, Hubbard campo médio (Sz={gs_sz})")
    data["hubbard_densidade_spin_fundamental"] = mf[gs_sz]["densidade_spin"]

    # 2. xTB: estados de spin candidatos = S de Lieb e S de Lieb + 1
    s_lieb = data["S_lieb"]
    mults = sorted({int(2 * s_lieb + 1), int(2 * s_lieb + 3)})
    data["xtb"] = {}
    best_geom = None
    for mult in mults:
        r, pos = xtb_optimize(symbols, coords, mult, freq=False)
        data["xtb"][str(mult)] = r
        write_xyz(out / f"geometria_xtb_mult{mult}.xyz", symbols, pos, f"{data['formula']} GFN2-xTB mult={mult}")
        print(f"  xTB mult={mult}: E={r['energia_eV']:.4f} eV, desvio do plano {r['desvio_do_plano_A']:.3f} Å", flush=True)
        if best_geom is None or r["energia_eV"] < best_geom[0]:
            best_geom = (r["energia_eV"], mult, pos)
    if not args.sem_freq:
        r, pos = xtb_optimize(symbols, best_geom[2], best_geom[1], fmax=0.01, freq=True)
        data["xtb_frequencias"] = {"multiplicidade": best_geom[1]} | r["frequencias"]
        print(f"  xTB frequências (mult={best_geom[1]}): imaginárias={r['frequencias']['imaginarias_cm1']}, "
              f"menor={r['frequencias']['menor_real_cm1']:.1f} cm⁻¹", flush=True)
        best_geom = (best_geom[0], best_geom[1], pos)

    # 3. DFT spin-polarizada na geometria xTB mais estável
    if not args.sem_dft:
        data["dft"] = {}
        base_2s = int(round(2 * s_lieb))
        runs = [(base_2s, False), (base_2s + 2, False)]
        if eta > abs(n_a - n_b):  # modos ocultos: o estado de spin baixo pode ter simetria quebrada
            runs.insert(1, (base_2s, True))
        for spin2s, bs in runs:
            r = dft_single_point(symbols, best_geom[2], spin2s, breaksym=bs)
            key = f"2S={spin2s}" + (" quebra-simetria" if bs else "")
            data["dft"][key] = r
            print(f"  DFT {key}: E={r['energia_Ha']:.6f} Ha, <S²>={r['S2']:.3f}, conv={r['convergiu']}", flush=True)
        e0 = min(r["energia_Ha"] for r in data["dft"].values())
        for r in data["dft"].values():
            r["energia_rel_meV"] = 27211.4 * (r["energia_Ha"] - e0)
        gs = min(data["dft"], key=lambda k: data["dft"][k]["energia_Ha"])
        data["dft_fundamental"] = gs
        spin_c = data["dft"][gs]["densidade_spin_atomos"][:n_c]
        spin_density_plot(out / "densidade_spin_dft.png", struct, spin_c,
                          f"{data['formula']} — densidade de spin, UB3LYP/6-31G ({gs})")
    data["tempo_s"] = time.time() - t0
    json.dump(data, open(out / "dados.json", "w"), indent=1, ensure_ascii=False)
    return data


REFERENCES = {
    "fenalenil": ([(0, 0), (1, 0), (0, 1)], "monorradical, S=1/2 (conhecido desde os anos 1950)"),
    "triangulena_3": ([(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (0, 2)],
                      "S=1 (tripleto) fundamental; sintetizada em superfície (Pavliček et al. 2017)"),
    "calice_de_clar": (None, "S=0 fundamental, tripleto 23 meV acima (Mishra et al., Nat. Nanotech. 2020)"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("candidato", nargs="?")
    ap.add_argument("--referencias", action="store_true")
    ap.add_argument("--sem-dft", action="store_true")
    ap.add_argument("--sem-freq", action="store_true")
    args = ap.parse_args()
    if args.referencias:
        for name, (cells, expected) in REFERENCES.items():
            if cells is None:  # o cálice de Clar é o C38H18 que as próprias civilizações acharam
                struct = json.load(open(ROOT / "data/candidatos/C38H18_eta2.json"))["estrutura"]
            else:
                struct = structure_from_cells(cells)
            run(name, struct, args, expected)
    if args.candidato:
        rec = json.load(open(args.candidato))
        struct = rec["estrutura"]
        # Confere a geometria reconstruída a partir dos anéis contra a exportada pelo C++.
        rebuilt = structure_from_cells([tuple(c) for c in struct["anel_coordenadas_axiais_q_r"]])
        assert len(rebuilt["carbonos"]) == len(struct["carbonos"]) and len(rebuilt["ligacoes"]) == len(struct["ligacoes"])
        run(pathlib.Path(args.candidato).stem, struct, args)


if __name__ == "__main__":
    main()
