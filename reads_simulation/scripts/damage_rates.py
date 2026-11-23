import gzip, sys
a, b = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 12
def seqs(p):
    with gzip.open(p, "rt") as f:
        for line in f:
            if line[0] != ">":
                yield line.strip()
c5 = [0]*n; t5 = [0]*n; g3 = [0]*n; a3 = [0]*n
for x, y in zip(seqs(a), seqs(b)):
    for i in range(min(n, len(x))):
        if x[i] == "C":
            c5[i] += 1
            if y[i] == "T": t5[i] += 1
        if x[-1-i] == "G":
            g3[i] += 1
            if y[-1-i] == "A": a3[i] += 1
print("pos\tC>T_5p\tG>A_3p")
for i in range(n):
    print("%d\t%.4f\t%.4f" % (i+1, t5[i]/max(c5[i],1), a3[i]/max(g3[i],1)))
