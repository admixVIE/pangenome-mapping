#!/bin/bash
#SBATCH --job-name=garg_sim
#SBATCH --array=0-8
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=08:00:00
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
scripts/01_simulate.sh $COND $DAM $CHUNK
