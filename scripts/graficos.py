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
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # ordem categorial fixa
CULTURE_COLOR = {"Empiristas": SERIES[0], "Incrementalistas": SERIES[1], "Evolucionistas": SERIES[2], "Isolados": SERIES[3]}
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
    for r in read(path):
        if r["melhor"] in ("nan", ""):
            continue
        best[r["cultura"]][int(r["ano"])].append(float(r["melhor"]))
        classes[r["cultura"]][int(r["ano"])].append(int(r["classes"]))
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, series, label in ((axes[0], best, ylabel), (axes[1], classes, "classes de fenômenos conhecidas")):
        for culture, color in CULTURE_COLOR.items():  # a cor segue a cultura, nunca a ordem do CSV
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


if __name__ == "__main__":
    emergence()
    research("brinquedo", "fração do ótimo global")
    research("quimica", "score (modos zero)")
    research("materiais", "score (planura + gap)")
    research("nuclear", "score (≈ Z do elemento mais pesado estável ≥ 60 s)")
    print("gráficos gerados em", RES)
