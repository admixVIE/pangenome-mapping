#!/bin/bash
# usage: 02_trim_merge.sh COND DAMAGE CHUNK
set -euo pipefail
source "$(dirname "$0")/config.sh"
ml load leeHom/1.2.18-GCC-14.3.0 >/dev/null 2>&1

COND=$1; DAM=$2; CHUNK=$3
IN=$SIM/$COND/$DAM/c$(printf '%03d' $CHUNK)
OUT=$TRIM/$COND/$DAM; mkdir -p $OUT
O=$OUT/c$(printf '%03d' $CHUNK)

[ -s ${O}.fq.gz ] && { echo "exists, skipping $O"; exit 0; }

leeHom --ancientdna \
  -fq1 ${IN}_s1.fq.gz -fq2 ${IN}_s2.fq.gz \
  -fqo $O \
  -t ${SLURM_CPUS_PER_TASK:-1} \
  --seed $(( SEED_BASE + CHUNK )) \
  --log ${O}.leehom.log
