#!/bin/bash
# usage: 01_simulate.sh COND DAMAGE CHUNK   (needs 01a_fragments.sh output; every stage is skipped if already finished)
set -euo pipefail
source "$(dirname "$0")/config.sh"
load_sim_env

COND=$1; DAM=$2; CHUNK=$3
SEED=$(( SEED_BASE + CHUNK ))
C=c$(printf '%03d' $CHUNK)
FRAG=$SIM/$COND/frags/$C.e.fa.gz
OUT=$SIM/$COND/$DAM; mkdir -p $OUT
P=$OUT/$C

[ -s ${P}_s2.fq.gz ] && { echo "exists, skipping $P"; exit 0; }
[ -s $FRAG ] || { echo "missing $FRAG, run 01a_fragments.sh first"; exit 1; }

if [ ! -s ${P}_d.fa.gz ]; then
  gzip -t $FRAG
  $GARG/src/deamSim --seed $(( SEED + 2 )) -damage ${DAMAGE[$DAM]} $FRAG | gzip > ${P}_d.tmp.fa.gz
  mv ${P}_d.tmp.fa.gz ${P}_d.fa.gz
fi

if [ ! -s ${P}_a.fa.gz ]; then
  $GARG/src/adptSim --seed $(( SEED + 3 )) -f $ADAPTF -s $ADAPTR -l $READLEN -artp ${P}_a.tmp.fa.gz ${P}_d.fa.gz
  mv ${P}_a.tmp.fa.gz ${P}_a.fa.gz
fi

timeout 6h $GARG/art_src_MountRainier/art_illumina -rs $(( SEED + 5 )) -ss $PLATFORM -amp -na -p \
  -i ${P}_a.fa.gz -l $READLEN -c 1 -qs 0 -qs2 0 -gzl 4 -o ${P}_t
mv ${P}_t1.fq.gz ${P}_s1.fq.gz
mv ${P}_t2.fq.gz ${P}_s2.fq.gz
rm -f ${P}_d.fa.gz ${P}_a.fa.gz
