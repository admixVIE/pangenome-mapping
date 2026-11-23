import gzip, itertools
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

N = 3000000
col = {"short": "#1f77b4", "mid": "#ff7f0e", "long": "#2ca02c"}
g = np.linspace(30, 240, 400)
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
for a, c in zip(ax, col):
    with gzip.open(f"reads/final/{c}/high/cov1/merged.fq.gz", "rt") as h:
        L = np.array([len(s) - 1 for s in itertools.islice(h, 1, 4 * N, 4)])
    a.plot(g, gaussian_kde(L, 0.1)(g), color=col[c], lw=2)
    a.axvline(np.mean(L), color="k", ls=":", lw=1)
    a.set_title(f"{c}: mean {L.mean():.1f}, median {np.median(L):.0f} bp (n={len(L):,})")
    a.set_xlabel("merged read length (bp)")
    a.set_xlim(30, 240)
ax[0].set_ylabel("density")
plt.tight_layout()
plt.savefig("figures/read_length_kde.png", dpi=160)
