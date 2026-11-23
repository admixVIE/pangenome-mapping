#!/bin/bash
#SBATCH --job-name=garg_artp
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the reads_simulation/ directory}"
source scripts/config.sh
I=${TID}
CHUNK=$(( I % NCHUNK ))
DAM=${DAMAGES[$(( (I / NCHUNK) % 3 ))]}
COND=${CONDS[$(( I / (NCHUNK * 3) ))]}
echo "== $COND $DAM chunk $CHUNK"
scripts/04_art_part.sh $COND $DAM $CHUNK $SLURM_ARRAY_TASK_ID
