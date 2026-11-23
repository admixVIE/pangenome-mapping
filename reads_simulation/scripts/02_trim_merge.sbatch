#!/bin/bash
#SBATCH --job-name=garg_trim
#SBATCH --array=0-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=24:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the reads_simulation/ directory}"
source scripts/config.sh
I=$SLURM_ARRAY_TASK_ID
CHUNK=$(( I % NCHUNK ))
DAM=${DAMAGES[$(( (I / NCHUNK) % 3 ))]}
COND=${CONDS[$(( I / (NCHUNK * 3) ))]}
echo "== $COND $DAM chunk $CHUNK"
scripts/02_trim_merge.sh $COND $DAM $CHUNK
