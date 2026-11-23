#!/bin/bash
#SBATCH --job-name=garg_frag
#SBATCH --array=0-2
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?submit from the reads_simulation/ directory}"
source scripts/config.sh
for CHUNK in $(seq 0 $(( NCHUNK - 1 ))); do scripts/01a_fragments.sh ${CONDS[$SLURM_ARRAY_TASK_ID]} $CHUNK; done
