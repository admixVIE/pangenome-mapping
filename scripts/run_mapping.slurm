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
#     bash scripts/map_to_pangenome.sh
# ==============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# --- Inputs (edit these) ---
S1="/path/to/reads_R1.fastq.gz"    # required
S2=""                               # optional second file for paired-end
THREADS="${SLURM_CPUS_PER_TASK:-8}"
KMER=29                             # vg k-mer length
WINDOW=11                           # vg window length

# --- Build argument list ---
args=(-s1 "$S1" -t "$THREADS" -k "$KMER" -w "$WINDOW")
[ -n "$S2" ] && args+=(-s2 "$S2")

echo "Running: map_to_pangenome.sh ${args[*]}"
bash "$SCRIPT_DIR/map_to_pangenome.sh" "${args[@]}"
