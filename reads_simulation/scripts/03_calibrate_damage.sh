#!/bin/bash
# usage: 03_calibrate_damage.sh NFRAG
set -euo pipefail
source "$(dirname "$0")/config.sh"
load_sim_env
stage_input
N=${1:-200000}
W=$ROOT/calibration; mkdir -p $W
$GARG/src/fragSim --seed 7 -n $N -m $MINSIZE -M $MAXSIZE -f $LEN/mid.freq $LOCAL_IN/endo/endo.fasta 2>/dev/null | gzip > $W/frag.fa.gz
while read -r name params; do
  $GARG/src/deamSim --seed 11 -damage $params $W/frag.fa.gz 2>/dev/null | gzip > $W/deam_$name.fa.gz
  echo "== $name  $params"
  python3 "$(dirname "$0")/damage_rates.py" $W/frag.fa.gz $W/deam_$name.fa.gz 6
done < <(for k in "${!DAMAGE[@]}"; do echo "$k ${DAMAGE[$k]}"; done | sort)
