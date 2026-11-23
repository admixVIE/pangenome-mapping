#!/bin/bash
#SBATCH --job-name=pangenome_map
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --output=map_%j.out
#SBATCH --error=map_%j.err

# ==============================================================================
# SLURM submission wrapper for map_to_pangenome.sh
# ------------------------------------------------------------------------------
# Edit the variables below (or pass them in via --export) and submit with:
#     sbatch scripts/run_mapping.slurm
#
# For usage / all options of the underlying script:
#     bash mapping/map_to_pangenome.sh
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Inputs (edit these) ---
SAMPLE="sample_name"                # required
RUNS_TSV="/path/to/runs.tsv"        # columns: sample, run, fq1, fq2 (fq2 empty for single-end)
THREADS="${SLURM_CPUS_PER_TASK:-8}"
KMER=15                             # mapping minimizer k (haplotype sampling is fixed at 29/11)
WINDOW=5                            # mapping minimizer w

echo "Running: map_to_pangenome.sh -S $SAMPLE --runs $RUNS_TSV -t $THREADS -k $KMER -w $WINDOW"
bash "$SCRIPT_DIR/map_to_pangenome.sh" -S "$SAMPLE" --runs "$RUNS_TSV" -t "$THREADS" -k "$KMER" -w "$WINDOW"
