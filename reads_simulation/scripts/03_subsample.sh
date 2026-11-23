#!/bin/bash
# usage: 03_subsample.sh COND DAMAGE
set -euo pipefail
source "$(dirname "$0")/config.sh"
ml load seqtk >/dev/null 2>&1

COND=$1; DAM=$2
SRC=$TRIM/$COND/$DAM
OUT=$FINAL/$COND/$DAM
mkdir -p $OUT/cov10

declare -A SUF=( [merged]="" [r1]="_r1" [r2]="_r2" )
for n in merged r1 r2; do
  [ -s $OUT/cov10/$n.fq.gz ] || cat $SRC/c[0-9][0-9][0-9]${SUF[$n]}.fq.gz > $OUT/cov10/$n.fq.gz
done

for COV in "${COVS[@]}"; do
  FRAC=$(python3 -c "print($COV/$COV_MAX)")
  mkdir -p $OUT/cov$COV
  pids=()
  for n in merged r1 r2; do
    [ -s $OUT/cov$COV/$n.fq.gz ] && continue
    ( seqtk sample -s $SUBSEED $OUT/cov10/$n.fq.gz $FRAC | gzip > $OUT/cov$COV/$n.fq.gz.tmp && mv $OUT/cov$COV/$n.fq.gz.tmp $OUT/cov$COV/$n.fq.gz ) &
    pids+=($!)
  done
  for p in "${pids[@]}"; do wait $p; done
  for n in merged r1 r2; do [ -s $OUT/cov$COV/$n.fq.gz ] || { echo "empty output cov$COV/$n"; exit 1; }; done
done
