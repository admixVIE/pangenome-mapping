#!/bin/bash
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p $INPUT/endo $INPUT/cont $INPUT/bact $LEN $ROOT/logs $SIM $TRIM

ln -sf $GENOME $INPUT/endo/endo.fasta
ln -sf $GENOME.fai $INPUT/endo/endo.fasta.fai

python3 - <<PY
import numpy as np
for name, ind in (("short","Ind1"),("mid","Ind2")):
    d = np.loadtxt("$IND_CDF_DIR/%s_cdf.txt" % ind)
    L, c = d[:,0].astype(int), d[:,1]
    f = np.diff(np.concatenate(([0.0], c)))
    f = f / f.sum()
    with open("$LEN/%s.freq" % name, "w") as o:
        for l, x in zip(L, f):
            o.write("%d\t%.10f\n" % (l, x))
    print(name, "mean", float((L*f).sum()))
PY
