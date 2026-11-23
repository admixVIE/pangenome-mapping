#!/bin/bash
# usage: 04_art_join.sh COND DAMAGE CHUNK
set -euo pipefail
source "$(dirname "$0")/config.sh"
COND=$1; DAM=$2; CHUNK=$3
C=c$(printf '%03d' $CHUNK)
OUT=$SIM/$COND/$DAM; PD=$OUT/parts
[ -s $OUT/${C}_s2.fq.gz ] && { echo "exists, skipping"; exit 0; }
NP=$(ls $PD | grep -c "^p[0-9][0-9]\.fa\.gz$")
for p in $(ls $PD | grep "^p[0-9][0-9]\.fa\.gz$" | sed 's/\.fa\.gz//'); do
  [ -s $PD/${p}_r1.fq.gz ] && [ -s $PD/${p}_r2.fq.gz ] || { echo "part $p incomplete"; exit 1; }
done
cat $PD/p[0-9][0-9]_r1.fq.gz > $OUT/${C}_s1.tmp.gz
cat $PD/p[0-9][0-9]_r2.fq.gz > $OUT/${C}_s2.tmp.gz
mv $OUT/${C}_s1.tmp.gz $OUT/${C}_s1.fq.gz
mv $OUT/${C}_s2.tmp.gz $OUT/${C}_s2.fq.gz
echo "joined $NP parts"
rm -rf $PD
rm -f $OUT/${C}_d.fa.gz $OUT/${C}_a.fa.gz $OUT/${C}_t1.fq.gz $OUT/${C}_t2.fq.gz
