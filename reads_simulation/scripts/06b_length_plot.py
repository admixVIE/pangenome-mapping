import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "pdf.fonttype": 42, "svg.fonttype": "none"})
col = {"short": "#0072B2", "mid": "#E69F00", "long": "#009E73"}
rows = []
fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4), sharex=True)
for a, c, tag in zip(ax, col, "abc"):
    d = np.loadtxt(f"figures/read_length_counts_{c}.tsv", dtype=np.int64)
    n = np.zeros(300, dtype=np.int64)
    n[d[:, 0]] = d[:, 1]
    L = np.arange(300)
    N = n.sum()
    p = n / N
    cum = np.cumsum(p)
    mean = (L * p).sum()
    sd = np.sqrt(((L - mean) ** 2 * p).sum())
    q = [int(np.searchsorted(cum, f)) for f in (0.05, 0.25, 0.5, 0.75, 0.95)]
    nz = np.nonzero(n)[0]
    rows.append((c, N, mean, sd, int(L[p.argmax()]), q, int(nz[0]), int(nz[-1])))
    a.fill_between(L, p, color=col[c], alpha=0.25, lw=0)
    a.plot(L, p, color=col[c], lw=1.4)
    a.set_title(f"Simulated {c} reads", fontsize=9)
    a.set_xlim(0, 250)
    a.set_ylim(bottom=0)
    a.yaxis.set_major_locator(MaxNLocator(4, steps=[1, 2, 5, 10]))
    a.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
    a.set_xlabel("Read length (bp)")
    a.text(-0.02, 1.12, tag, transform=a.transAxes, fontsize=11, fontweight="bold", va="bottom", ha="right")
    for s in ("top", "right"):
        a.spines[s].set_visible(False)
    a.tick_params(direction="out", length=3)
ax[0].set_ylabel("Density")
fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(f"figures/read_length_distribution.{ext}", dpi=600)
with open("figures/read_length_summary.tsv", "w") as o:
    print("condition\tn_merged_reads\tmean\tsd\tmode\tp5\tp25\tmedian\tp75\tp95\tmin\tmax", file=o)
    for c, N, m, sd, mo, q, lo, hi in rows:
        print(f"{c}\t{N}\t{m:.2f}\t{sd:.2f}\t{mo}\t" + "\t".join(map(str, q)) + f"\t{lo}\t{hi}", file=o)
print(open("figures/read_length_summary.tsv").read())
