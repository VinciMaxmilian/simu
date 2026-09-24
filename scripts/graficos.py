"""Gera os gráficos a partir dos CSVs em resultados/.

Uso: python3 scripts/graficos.py
"""
import csv
import pathlib
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
RES = ROOT / "resultados"
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]  # ordem categorial fixa
INK, MUTED, GRID = "#1f1f1e", "#6b6a64", "#e4e3dc"

plt.rcParams.update({
    "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": MUTED, "xtick.color": MUTED,
    "ytick.color": MUTED, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
    "axes.spines.top": False, "axes.spines.right": False, "axes.titlecolor": INK,
    "axes.titlesize": 11, "legend.frameon": False,
})


def read(path):
    with open(path) as fh:
        return list(csv.DictReader(fh))


def emergence():
    path = RES / "emergencia/amostras.csv"
    if not path.exists():
        return
    data = defaultdict(list)
    for r in read(path):
        data[r["cenario"]].append(float(r["tempo_apos_formacao_gyr"]))
    fig, ax = plt.subplots(figsize=(8, 4.2))
    for i, (name, vals) in enumerate(data.items()):
        ax.hist(vals, bins=80, range=(0, 20), density=True, histtype="step", linewidth=2, color=SERIES[i],
                label=name)
    ax.axvline(4.54, color=INK, linewidth=1, linestyle="--")
    ax.text(4.64, ax.get_ylim()[1] * 0.92, "Terra (4,54 Gyr)", color=INK)
    ax.set_xlabel("bilhões de anos desde a formação do planeta até o início da ciência")
    ax.set_ylabel("densidade (planetas que chegaram lá)")
    ax.set_title("Quando surge uma civilização científica?", loc="left")
    ax.legend()
    fig.tight_layout()
    fig.savefig(RES / "emergencia/grafico.png", dpi=150)
    plt.close(fig)


def research(lab, ylabel):
    path = RES / lab / "linha_do_tempo.csv"
    if not path.exists():
        return
    best, classes = defaultdict(lambda: defaultdict(list)), defaultdict(lambda: defaultdict(list))
    names = {}  # índice da civilização (ordem fixa da lista) → nome; a cor segue o índice
    for r in read(path):
        names[int(r["civilizacao"])] = r["cultura"]
        if r["melhor"] in ("nan", ""):
            continue
        best[r["cultura"]][int(r["ano"])].append(float(r["melhor"]))
        classes[r["cultura"]][int(r["ano"])].append(int(r["classes"]))
    culture_color = {names[i]: SERIES[i % len(SERIES)] for i in sorted(names)}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, series, label in ((axes[0], best, ylabel), (axes[1], classes, "classes de fenômenos conhecidas")):
        for culture, color in culture_color.items():  # a cor segue a civilização, nunca a ordem do CSV
            by_year = series.get(culture)
            if not by_year:
                continue
            years = sorted(by_year)
            means = [sum(by_year[y]) / len(by_year[y]) for y in years]
            ax.plot(years, means, color=color, linewidth=2, label=culture)
        ax.set_xlabel("anos desde o início da ciência")
        ax.set_ylabel(label)
    axes[0].set_title("Melhor resultado (média entre universos)", loc="left")
    axes[1].set_title("Conhecimento acumulado", loc="left")
    axes[1].legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(RES / lab / "grafico.png", dpi=150)
    plt.close(fig)


def biochemistry():
    path = RES / "bioquimica/amostras.csv"
    if not path.exists():
        return
    rows = [r for r in read(path) if r["cenario"] == "Terra" and r["catalise"] == "adaptada"]
    by = defaultdict(list)
    for r in rows:
        by[r["bioquimica"]].append(float(r["tempo_apos_formacao_gyr"]))
    order = sorted(by, key=lambda k: sorted(by[k])[len(by[k]) // 2])
    fig, ax = plt.subplots(figsize=(8, 0.45 * len(order) + 1.4))
    for i, name in enumerate(order):
        v = sorted(by[name])
        q = lambda f: v[int(f * (len(v) - 1))]
        ax.plot([q(0.1), q(0.9)], [i, i], color=SERIES[0], linewidth=2, solid_capstyle="round")
        ax.plot([q(0.5)], [i], "o", color=SERIES[0], markersize=8)
        ax.text(q(0.9) + 0.2, i, f"{q(0.5):.1f}", va="center", color=INK, fontsize=8)
    ax.axvline(4.54, color=INK, linewidth=1, linestyle="--")
    ax.set_yticks(range(len(order)), order)
    ax.set_xlabel("bilhões de anos da formação do planeta até a ciência (mediana e faixa 10–90%)")
    ax.set_title("Cenário Terra, catálise adaptada: quando cada bioquímica chega à ciência", loc="left")
    ax.grid(axis="y", visible=False)
    fig.tight_layout()
    fig.savefig(RES / "bioquimica/grafico.png", dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    emergence()
    research("brinquedo", "fração do ótimo global")
    research("quimica", "score (modos zero)")
    research("materiais", "score (planura + gap)")
    research("nuclear", "score (≈ Z do elemento mais pesado estável ≥ 60 s)")
    biochemistry()
    for lab, label in (("brinquedo", "fração do ótimo global"), ("quimica", "score (modos zero)"),
                       ("materiais", "score (planura + gap)")):
        research(lab + "_bioquimicas", label)
    print("gráficos gerados em", RES)
