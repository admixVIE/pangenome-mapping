import gzip, itertools, subprocess, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

R = "reads/final/{c}/high/cov1/{f}.fq.gz"
CONDS = ["short", "mid", "long"]
N = int(sys.argv[1]) if len(sys.argv) > 1 else 3000000
LO, HI = 30, 291
x = np.arange(LO, HI + 1)

def target(c):
    if c == "long":
        p = np.exp(-(np.log(x) - 4.6) ** 2 / (2 * 0.4 ** 2)) / x
    else:
        p = np.loadtxt(f"length/{c}.freq")[:, 1]
    return p / p.sum()

def read_stats(c):
    L, F = [], []
    with gzip.open(R.format(c=c, f="merged"), "rt") as h:
        for name, seq in ((l, s) for l, s in zip(itertools.islice(h, 0, 4 * N, 4), itertools.islice(h, 0, None, 4))):
            pass
    return L, F

def load(c):
    L, F = [], []
    with gzip.open(R.format(c=c, f="merged"), "rt") as h:
        for i, line in enumerate(itertools.islice(h, 0, 4 * N)):
            k = i % 4
            if k == 0:
                parts = line.split(":")
                F.append(int(parts[3]) - int(parts[2]))
            elif k == 1:
                L.append(len(line) - 1)
    return np.array(L), np.array(F)

def count(c, f):
    out = subprocess.run(f"zcat {R.format(c=c, f=f)} | wc -l", shell=True, capture_output=True, text=True).stdout
    return int(out) // 4

col = {"short": "#1f77b4", "mid": "#ff7f0e", "long": "#2ca02c"}
fig, ax = plt.subplots(1, 3, figsize=(15, 4), sharey=False)
rows = []
for a, c in zip(ax, CONDS):
    L, F = load(c)
    t = target(c)
    tm = (x * t).sum()
    tmed = x[np.searchsorted(np.cumsum(t), 0.5)]
    nm, nu = count(c, "merged"), count(c, "r1")
    g = np.linspace(LO, HI, 400)
    a.plot(g, gaussian_kde(L, 0.1)(g), color=col[c], lw=2, label="merged reads")
    a.plot(g, gaussian_kde(F[F <= HI], 0.1)(g), color="k", lw=1.2, ls="--", label="fragments (read names)")
    a.plot(x, t / (x[1] - x[0]), color="r", lw=1.2, ls=":", label="target")
    a.set_title(f"{c} (n={len(L):,})")
    a.set_xlabel("length (bp)")
    a.set_xlim(LO, HI)
    a.legend(frameon=False)
    ks = np.max(np.abs(np.cumsum(np.bincount(L, minlength=HI + 1)[LO:HI + 1]) / len(L) - np.cumsum(t)))
    rows.append((c, tm, tmed, F.mean(), np.median(F), L.mean(), np.median(L), L.min(), L.max(), nm, nu, nm / (nm + nu), ks))
ax[0].set_ylabel("density")
plt.tight_layout()
plt.savefig("figures/read_length_kde.png", dpi=160)
hdr = "cond\ttarget_mean\ttarget_median\tfrag_mean\tfrag_median\tmerged_mean\tmerged_median\tmin\tmax\tn_merged\tn_unmerged_pairs\tfrac_merged\tKS_vs_target"
with open("figures/read_length_table.tsv", "w") as o:
    print(hdr, file=o)
    for r in rows:
        print("\t".join(r[0:1] + tuple(f"{v:.2f}" if isinstance(v, float) else str(v) for v in r[1:])), file=o)
print(open("figures/read_length_table.tsv").read())
