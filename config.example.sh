# ==============================================================================
# config.example.sh  —  TEMPLATE.  Copy it to config.sh and edit for YOUR system:
#
#       cp config.example.sh config.sh
#       nano config.sh
#
# ------------------------------------------------------------------------------
# WHY THIS FILE EXISTS
# The pipeline scripts contain NO machine-specific paths. Everything that is
# different from one computer to another — where your reference files live and how the required tools are made available — lives
# here in config.sh. map_to_pangenome.sh loads config.sh automatically before it
# runs, so anything you set below overrides the built-in defaults.
#
# config.sh is git-ignored, so YOUR local paths are never committed. Everyone who
# uses the repo keeps their own private config.sh; this template stays shared.
#
# ON ANY MACHINE OTHER THAN THE ORIGINAL LISC CLUSTER, CREATING config.sh IS
# REQUIRED. Without it the scripts fall back to LISC-specific module and conda
# names that will not exist on your system.
#
# You must set TWO things correctly (everything else usually stays as-is):
#   (a) the paths to your reference files, and
#   (b) load_environment(), which makes the tools available on your system.
# ==============================================================================


# ------------------------------------------------------------------------------
# (a) PATHS
# ------------------------------------------------------------------------------

# Directory where references/download_references.sh downloaded the big files.
# Edit this to the absolute path of your references/ directory.
export PANGENOME_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/references"

# Pangenome graph fileset. Leave as-is for the HPRC v1.1 graph shipped by
# download_references.sh. Only change these if you use a different graph.
export GRAPH_NAME="hprc-v1.1-mc-grch38"          # expects $GRAPH_NAME.gbz
export HAPL_NAME="hprc-v1.1-mc-grch38.hapl"
export PATHS_NAME="hprc-v1.1-mc-grch38.paths_list.txt"
export SEQ_PREFIX="GRCh38#0#"                     # contig-name prefix to strip

# Linear reference genome. It must have a .fai and a .dict file next to it
# (download_references.sh builds both for you).
export REF_FA="$PANGENOME_ROOT/hg38.fa"

# Mapping minimizer (used by vg minimizer / giraffe). Optional: if unset, vg
# giraffe's own defaults apply (k=29, w=11). The -k / -w command-line flags
# override these.
# export MAP_KMER=29
# export MAP_WINDOW=11

# Number of candidate haplotypes for personalized graph sampling (default 4).
# Override per run with --num-haplotypes.
# export NUM_HAPLOTYPES=4


# ------------------------------------------------------------------------------
# (b) ENVIRONMENT
# ------------------------------------------------------------------------------
# load_environment() is run for you, once, before the pipeline starts. Put here
# whatever commands YOUR system needs so that these tools end up on the PATH:
#
#       vg   kmc   samtools   python3
#
# Pick the block that matches your setup and delete the others.

load_environment() {

    # --- Option 1: an environment-modules / Lmod cluster ---
    # (module names are examples — use whatever `module avail` shows on your cluster)
    ml load KMC vg SAMtools

    # --- Option 2: a single conda environment with everything installed ---
    # conda activate pangenome-mapping

    # --- Option 3: manual installs on a plain server ---
    # export PATH="/opt/vg/bin:/opt/kmc/bin:$PATH"

    :  # (no-op so the function is never empty; remove once you add commands above)
}
