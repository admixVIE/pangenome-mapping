# ==============================================================================
# config.example.sh  ->  copy to config.sh and edit for YOUR system.
# ------------------------------------------------------------------------------
# map_to_pangenome.sh sources config.sh (if present) before running. Everything
# here is optional: with no config.sh the script falls back to the LISC cluster
# defaults it was developed on. On any other machine, set at least PANGENOME_ROOT
# and load_environment().
#
# config.sh is git-ignored, so your local paths never get committed.
# ==============================================================================

# --- Where download_references.sh put the reference files ---
export PANGENOME_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/references"

# --- Pangenome graph fileset (usually leave as-is for the HPRC v1.1 graph) ---
export GRAPH_NAME="hprc-v1.1-mc-grch38"          # expects $GRAPH_NAME.gbz
export HAPL_NAME="hprc-v1.1-mc-grch38.hapl"
export PATHS_NAME="hprc-v1.1-mc-grch38.paths_list.txt"
export SEQ_PREFIX="GRCh38#0#"                     # prefix stripped from contig names

# --- Reference genome (must match the graph's reference paths) ---
export REF_FA="$PANGENOME_ROOT/hg38.fa"           # needs .fai and .dict alongside

# --- Tools ---
export ABRA2_JAR="$PANGENOME_ROOT/abra2-2.23.jar"
# GATK 3.8 RealignerTargetCreator needs Java 8:
export GATK3_JAR="/path/to/GenomeAnalysisTK.jar"  # GATK 3.8 jar
export JAVA8="java"                               # a Java 8 executable

# --- Make the pipeline's tools available (edit for your cluster/conda) ---
# Required on PATH: kmc, vg, samtools, bedtools, bamleftalign (freebayes), java.
load_environment() {
    # Example for an environment-modules cluster:
    ml load KMC vg Conda freebayes SAMtools BEDTools
    conda activate GATK-3.8.0
    conda activate freebayes
    # Or, for a single conda env:
    # conda activate pangenome-mapping
}
