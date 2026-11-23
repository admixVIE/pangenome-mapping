#!/usr/bin/env python3
# summary/mapdamage/{low,mid,high}/{short,mid,long}/{5pCtoT,3pGtoA}_freq.txt (mapDamage on linear-arm 10x BAMs) -> summary/mapdamage/mapdamage_combined.{png,pdf}
# lines = damage level (one sequential blue hue, light -> dark), value = mean over the three read-length classes (they differ by < 0.01).
# run:  ml purge; ml load matplotlib/3.7.2-gfbf-2023a; python3 plot_mapdamage_combined.py
import os, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
B = os.environ.get("BENCH_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))); D = f"{B}/summary/mapdamage"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e3e2de", "#fcfcfb"
DAM = [("low", "#6da7ec", "o"), ("mid", "#2a78d6", "s"), ("high", "#184f95", "^")]
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF})
fig, axs = plt.subplots(1, 2, figsize=(10, 4.3), sharey=True)
for ax, (f, col, ttl) in zip(axs, [("5pCtoT_freq.txt", "5pC>T", "5' end: C→T"), ("3pGtoA_freq.txt", "3pG>A", "3' end: G→A")]):
    ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
    for sp in ("top", "right"): ax.spines[sp].set_visible(False)
    for dm, c, mk in DAM:
        t = [pd.read_csv(f"{D}/{dm}/{L}/{f}", sep="\t").set_index("pos")[col] for L in ("short", "mid", "long")]
        y = sum(t) / 3
        ax.plot(y.index, y.values, color=c, lw=2.0, marker=mk, ms=4.5, mec=SURF, mew=0.7, label=f"{dm} damage")
        print(ttl, dm, "pos1 =", round(y.iloc[0], 3))
    ax.set_title(ttl, loc="left", fontsize=11.5); ax.set_xlabel("position from the read end (bp)")
axs[0].set_ylabel("substitution frequency"); axs[0].legend(frameon=False, loc="upper right")
fig.suptitle("mapDamage profiles of the three simulated damage levels", x=0.01, ha="left", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for ext in ("png", "pdf"): fig.savefig(f"{D}/mapdamage_combined.{ext}", dpi=200 if ext == "png" else None, facecolor=SURF)
print("saved", f"{D}/mapdamage_combined.png")
