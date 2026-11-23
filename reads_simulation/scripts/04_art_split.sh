#!/bin/bash
# usage: 04_art_split.sh COND DAMAGE CHUNK   (optional env NREC overrides the record count)
set -euo pipefail
source "$(dirname "$0")/config.sh"
COND=$1; DAM=$2; CHUNK=$3
C=c$(printf '%03d' $CHUNK)
OUT=$SIM/$COND/$DAM; PD=$OUT/parts
NTOTAL=$(python3 -c "print(round($COV_MAX*$GENOME_BP/${MEANLEN[$COND]}))")
N=${NREC:-$(( (NTOTAL + NCHUNK - 1) / NCHUNK ))}
L=$(( 2 * ( (N + ART_PARTS - 1) / ART_PARTS ) ))
[ -s $OUT/${C}_a.fa.gz ] || { echo "missing $OUT/${C}_a.fa.gz"; exit 1; }
[ -e $PD/.split_done ] && { echo "already split"; exit 0; }
rm -rf $PD; mkdir -p $PD
zcat $OUT/${C}_a.fa.gz | split -d -a 2 -l $L --filter='gzip -1 > "$FILE.fa.gz"' - $PD/p
ls $PD | grep -c "fa.gz$"
touch $PD/.split_done
