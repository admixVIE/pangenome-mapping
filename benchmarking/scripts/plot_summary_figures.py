#!/usr/bin/env python3
# summary/linear_vs_pangenome_summary.tsv -> summary/figures/{precision,recall,F1}.{png,pdf}
# one figure per metric; rows = SNP / INDEL, columns = read length; x = coverage (DP filter in tick label);
# colour = mapping arm (light blue pangenome, pink linear; validated, contrast WARN -> relief via legend + line styles), line style + marker = damage.
# run:  ml purge; ml load matplotlib/3.7.2-gfbf-2023a; python3 plot_summary_figures.py
import os, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = os.environ.get("BENCH_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
df = pd.read_csv(f"{B}/summary/linear_vs_pangenome_summary.tsv", sep="\t")
out = f"{B}/summary/figures"; os.makedirs(out, exist_ok=True)
BLUE, PINK, INK, INK2, GRID, SURF = "#6da7ec", "#e87ba4", "#0b0b0b", "#52514e", "#e3e2de", "#fcfcfb"
ARMS = [("pangenome", BLUE), ("linear", PINK)]
DAM = [("low", "-", "o"), ("mid", "--", "s"), ("high", ":", "^")]
COV = [("1x", "DP>=1"), ("5x", "DP>=5"), ("7x", "DP>=7"), ("10x", "DP>=7")]
XL = [f"{c}\n{d.replace('>=', '≥')}" for c, d in COV]
LEN = ["short", "mid", "long"]
METRICS = [("precision", "Precision"), ("recall", "Recall"), ("F1", "F1")]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF})
for m, mname in METRICS:
    fig, axs = plt.subplots(2, 3, figsize=(11.5, 7.2), sharex=True)
    for r, t in enumerate(["SNP", "INDEL"]):
        lo = df[df["type"] == t][[f"{a}_{m}" for a, _ in ARMS]].min().min(); hi = df[df["type"] == t][[f"{a}_{m}" for a, _ in ARMS]].max().max()
        pad = (hi - lo) * 0.08
        for c, ln in enumerate(LEN):
            ax = axs[r, c]; ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
            for sp in ("top", "right"): ax.spines[sp].set_visible(False)
            for arm, col in ARMS:
                for dm, ls, mk in DAM:
                    ys = []
                    for cv, dp in COV:
                        s = df[(df["type"] == t) & (df["length"] == ln) & (df["damage"] == dm) & (df["cov"] == cv) & (df["DP_filter"] == dp)]
                        ys.append(float(s[f"{arm}_{m}"].iloc[0]) if len(s) else float("nan"))
                    ax.plot(range(len(COV)), ys, color=col, ls=ls, marker=mk, ms=5, lw=2.0, mec=SURF, mew=0.8)
            ax.set_ylim(max(0, lo - pad), min(1, hi + pad)); ax.set_xticks(range(len(COV))); ax.set_xticklabels(XL, fontsize=9)
            ax.set_title(f"{t} - {ln} reads", fontsize=11, color=INK, loc="left")
            if c == 0: ax.set_ylabel(mname)
    h = [Line2D([], [], color=col, lw=2.2, label=a) for a, col in ARMS] + [Line2D([], [], color=INK2, ls=ls, marker=mk, ms=5, lw=1.8, label=f"{dm} damage") for dm, ls, mk in DAM]
    fig.legend(handles=h, loc="upper center", ncol=5, frameon=False, bbox_to_anchor=(0.5, 0.965), fontsize=10)
    fig.suptitle(f"{mname}: pangenome vs linear mapping", fontsize=12.5, x=0.01, ha="left", y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.925))
    for ext in ("png", "pdf"): fig.savefig(f"{out}/{m}.{ext}", dpi=200 if ext == "png" else None, facecolor=SURF)
    plt.close(fig); print("saved", f"{out}/{m}.png")
