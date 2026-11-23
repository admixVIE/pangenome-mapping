#!/usr/bin/env python3
# summary/read_lengths/hist_{short,mid,long}.tsv (length<TAB>count; first 20M reads of each low-damage 5x set) -> summary/read_lengths/read_length_kde.{png,pdf} + read_length_stats.tsv
# run:  ml purge; ml load matplotlib/3.7.2-gfbf-2023a; python3 plot_read_length_kde.py
import os, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
B = os.environ.get("BENCH_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__))); D = f"{B}/summary/read_lengths"
INK, INK2, GRID, SURF = "#0b0b0b", "#52514e", "#e3e2de", "#fcfcfb"
COL = {"short": "#2a78d6", "mid": "#eb6834", "long": "#1baf7a"}     # categorical slots 1-3 of the validated palette
plt.rcParams.update({"font.size": 10, "axes.edgecolor": INK2, "axes.labelcolor": INK, "text.color": INK, "xtick.color": INK2, "ytick.color": INK2, "figure.facecolor": SURF, "axes.facecolor": SURF})
H = {L: pd.read_csv(f"{D}/hist_{L}.tsv", sep="\t", header=None, names=["len", "n"]) for L in COL}
stats = []
for L, h in H.items():
    n = h["n"].sum(); cum = h["n"].cumsum() / n; q = lambda p: int(h["len"][np.searchsorted(cum.values, p)])
    mean = float((h["len"] * h["n"]).sum() / n)
    stats.append(dict(set=L, reads_counted=int(n), min=int(h["len"].min()), q05=q(.05), median=q(.5), mean=round(mean, 1), q95=q(.95), max=int(h["len"].max())))
st = pd.DataFrame(stats); st.to_csv(f"{D}/read_length_stats.tsv", sep="\t", index=False); print(st.to_string(index=False))
fig, ax = plt.subplots(figsize=(8.5, 4.6)); ax.grid(axis="y", color=GRID, lw=0.8); ax.set_axisbelow(True)
for sp in ("top", "right"): ax.spines[sp].set_visible(False)
xmax = max(s["q95"] for s in stats) * 1.25
x = np.linspace(0, xmax, 1200)
for L, h in H.items():
    kde = gaussian_kde(h["len"].values.astype(float), weights=h["n"].values.astype(float), bw_method=0.05)   # bw = 5% of the std
    y = kde(x); ax.plot(x, y, color=COL[L], lw=2.2, label=f"{L}  (median {st.loc[st.set == L, 'median'].iloc[0]} bp)")
    ax.fill_between(x, y, color=COL[L], alpha=0.12, lw=0)
ax.set_xlim(0, xmax); ax.set_ylim(bottom=0); ax.set_xlabel("read length (bp, merged reads)"); ax.set_ylabel("density")
ax.legend(frameon=False, loc="upper right"); ax.set_title("Read-length distribution of the simulated reads", loc="left", fontsize=12)
fig.tight_layout()
for ext in ("png", "pdf"): fig.savefig(f"{D}/read_length_kde.{ext}", dpi=200 if ext == "png" else None, facecolor=SURF)
print("saved", f"{D}/read_length_kde.png")
