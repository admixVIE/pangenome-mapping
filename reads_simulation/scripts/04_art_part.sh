#!/bin/bash
# usage: 04_art_part.sh COND DAMAGE CHUNK PART
set -euo pipefail
source "$(dirname "$0")/config.sh"
load_sim_env
COND=$1; DAM=$2; CHUNK=$3; PART=$4
SEED=$(( SEED_BASE + CHUNK ))
PD=$SIM/$COND/$DAM/parts
q=$(printf '%02d' $PART)
[ -s $PD/p${q}_r2.fq.gz ] && { echo "part $q done"; exit 0; }
[ -s $PD/p$q.fa.gz ] || { echo "no part $q"; exit 0; }

for ATT in $(seq 0 $(( ART_ATTEMPTS - 1 ))); do
  rm -f $PD/p${q}_o1.fq.gz $PD/p${q}_o2.fq.gz
  S=$(( SEED + 5 + 1000 * (PART + 1) + 7 * ATT ))
  if timeout $ART_PART_TIMEOUT $GARG/art_src_MountRainier/art_illumina -rs $S -ss $PLATFORM -amp -na -p \
       -i $PD/p$q.fa.gz -l $READLEN -c 1 -qs 0 -qs2 0 -gzl 4 -o $PD/p${q}_o > $PD/p$q.art.log 2>&1; then
    mv $PD/p${q}_o1.fq.gz $PD/p${q}_r1.fq.gz
    mv $PD/p${q}_o2.fq.gz $PD/p${q}_r2.fq.gz
    echo "part $q ok on attempt $ATT seed $S"
    exit 0
  fi
  echo "part $q attempt $ATT (seed $S) failed or timed out"
done
exit 1
