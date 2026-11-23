#!/bin/bash
#SBATCH --job-name=garg_cal
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
cd "${SLURM_SUBMIT_DIR:?submit from the reads_simulation/ directory}"
scripts/03_calibrate_damage.sh 500000
