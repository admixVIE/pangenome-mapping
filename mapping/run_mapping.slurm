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
#     sbatch mapping/run_mapping.slurm
#
# For usage / all options of the underlying script:
#     bash mapping/map_to_pangenome.sh
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Inputs (edit these) ---
SAMPLE="sample_name"                # required
RUNS_TSV="/path/to/runs.tsv"        # columns: sample, run, fq1, fq2 (fq2 empty for single-end)
THREADS="${SLURM_CPUS_PER_TASK:-8}"
KMER=""                             # mapping minimizer k; empty = config.sh MAP_KMER or giraffe default (29)
WINDOW=""                           # mapping minimizer w; empty = config.sh MAP_WINDOW or giraffe default (11)

args=(-S "$SAMPLE" --runs "$RUNS_TSV" -t "$THREADS")
[ -n "$KMER" ]   && args+=(-k "$KMER")
[ -n "$WINDOW" ] && args+=(-w "$WINDOW")

echo "Running: map_to_pangenome.sh ${args[*]}"
bash "$SCRIPT_DIR/map_to_pangenome.sh" "${args[@]}"
