#!/bin/bash
# usage: 01a_fragments.sh COND CHUNK [NFRAG_OVERRIDE]
set -euo pipefail
source "$(dirname "$0")/config.sh"
load_sim_env
stage_input

COND=$1; CHUNK=$2
NTOTAL=$(python3 -c "print(round($COV_MAX*$GENOME_BP/${MEANLEN[$COND]}))")
N=${3:-$(( (NTOTAL + NCHUNK - 1) / NCHUNK ))}
SEED=$(( SEED_BASE + CHUNK ))
D=$SIM/$COND/frags; mkdir -p $D
F=$D/c$(printf '%03d' $CHUNK).e.fa.gz
[ -s $F ] && { echo "exists, skipping $F"; exit 0; }

$GARG/src/fragSim --seed $(( SEED + 1 )) -tag e -n $N -m $MINSIZE -M $MAXSIZE ${LENARGS[$COND]} $LOCAL_IN/endo/endo.fasta | gzip > $F.tmp.gz
mv $F.tmp.gz $F
