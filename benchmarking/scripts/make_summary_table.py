#!/usr/bin/env python3
# Rerunnable: postproc/results/summary_happy.tsv -> summary/linear_vs_pangenome_summary.{tsv,md}
# one row per (type, length, damage, coverage, callset); pangenome and linear side by side. Missing arm/set = NA.
import csv, os
B = os.environ.get("BENCH_ROOT") or os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = f"{B}/postproc/results/summary_happy.tsv"; out = f"{B}/summary/linear_vs_pangenome_summary"
d = {}
for x in csv.DictReader(open(src), delimiter="\t"):
    if x["callset"] == "confidence" and x["cov"] == "1": continue   # 1x at DP>=5 dropped (recall ~0, uninformative); 1x DP>=1 (confidence_DP1) is kept
    d[(x["Type"], x["set"], x["callset"], x["arm"])] = x
dp = {"1": "DP>=5", "5": "DP>=5", "7": "DP>=7", "10": "DP>=7"}
def key(k):
    s = k[1]; ln, dm, cv = s.split("_"); cv = int(cv[3:])
    return (k[0] != "SNP", ["short", "mid", "long"].index(ln), cv, k[2] != "confidence", ["low", "mid", "high"].index(dm))
rows = sorted({(t, s, c) for t, s, c, a in d}, key=key)
M = ["METRIC.Precision", "METRIC.Recall", "METRIC.F1_Score"]
hdr = ["type", "length", "damage", "cov", "DP_filter"]
for a in ("pangenome", "linear"): hdr += [f"{a}_precision", f"{a}_recall", f"{a}_F1"]
hdr += ["F1_diff_pangenome_minus_linear"]
def f(v): return "NA" if v in (None, "") else f"{float(v):.3f}"
tab = []
for t, s, c in rows:
    ln, dm, cv = s.split("_"); cv = cv[3:]
    r = [t, ln, dm, cv + "x", "DP>=1" if c == "confidence_DP1" else dp[cv]]
    v = {}
    for a in ("pangenome", "linear"):
        x = d.get((t, s, c, a)); v[a] = float(x["METRIC.F1_Score"]) if x else None
        r += [f(x[m]) if x else "NA" for m in M]
    r.append("NA" if None in v.values() else f"{v['pangenome']-v['linear']:+.3f}")
    tab.append(r)
with open(out + ".tsv", "w") as o:
    o.write("\t".join(hdr) + "\n"); [o.write("\t".join(r) + "\n") for r in tab]
with open(out + ".md", "w") as o:
    o.write("| " + " | ".join(hdr) + " |\n|" + "---|" * len(hdr) + "\n"); [o.write("| " + " | ".join(r) + " |\n") for r in tab]
print(len(tab), "rows ->", out + ".tsv /.md")
