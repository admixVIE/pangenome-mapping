# benchmarking/config.example.sh -- copy to benchmarking/config.sh (git-ignored) and edit.
# Sourced by happy_array.sbatch / happy_dp1_array.sbatch.

# Working directory holding the per-arm results (postproc/results/<arm>/<set>/ ...).
# The python scripts read BENCH_ROOT too (export it before running them).
export BENCH_ROOT=/path/to/benchmark_workdir
B=$BENCH_ROOT

# This repo's benchmarking/lib (lib_harmonize.sh, hetfilt awk, rtg-wrapped-mem).
BP=/path/to/pangenome-mapping/benchmarking/lib

# Truth set directory: GRCh38_HG2-T2TQ100-V1.1_smvar.{vcf.gz,benchmark.bed},
# chrom_map.tsv, primary_contigs.bed, hg38.sdf, stratifications/.
T_=/path/to/truth

REF_FA=/path/to/hg38.fa
HAPPY_SIF=/path/to/happy.sif          # hap.py apptainer image
BIND_PATHS=/path/to:/path/to          # apptainer --bind (host:container)

# rtg-tools (used by hap.py's vcfeval engine, via lib/rtg-wrapped-mem)
export RTG_DIR=/path/to/rtg-tools
export RTG_JAVA_BIN=java
