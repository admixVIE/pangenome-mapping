#!/bin/bash
# hap.py (genotype-level, as in with_wave_03 / GIAB README) on confidence.vcf.gz of one finished set; rtg heap capped (RTG_MEM).
#   list line N of $LIST = "<arm> <set>"
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=08:00:00
set -uo pipefail
# Paths come from benchmarking/config.sh (copy config.example.sh); submit from benchmarking/ or set BENCH_CONFIG.
source "${BENCH_CONFIG:-${SLURM_SUBMIT_DIR:?submit from benchmarking/}/config.sh}"
read ARM SET < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ${LIST:?LIST not set})
COV=${SET##*cov}; W=$B/postproc/results/$ARM/$SET; HD=$W/happy/confidence
[ -e $W/.done_sets ] && [ -s $W/confidence.vcf.gz ] || { echo "not finished: $ARM $SET"; exit 1; }
[ -e $HD/.done ] && { echo "already scored $ARM $SET"; exit 0; }
ml load BCFtools SAMtools >/dev/null 2>&1; source $BP/lib_harmonize.sh
mkdir -p $HD; rm -f $HD/happy.*; echo "== hap.py $ARM $SET $(date)"
harmonize_query $W/confidence.vcf.gz $HD/harm $T_/chrom_map.tsv $T_/primary_contigs.bed || exit 1
RTGMEM=64g apptainer exec --bind "$BIND_PATHS" --env RTG_TALKBACK=false $HAPPY_SIF \
  hap.py --engine vcfeval --engine-vcfeval-path $BP/rtg-wrapped-mem --engine-vcfeval-template $T_/hg38.sdf \
    --gender male -r $REF_FA -f $T_/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed \
    --stratification $T_/stratifications/GRCh38@all/GRCh38-all-stratifications.abspath.tsv \
    -o $HD/happy --threads $SLURM_CPUS_PER_TASK $T_/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz $HD/harm.vcf.gz > $HD/happy.log 2>&1
if [ -s $HD/happy.summary.csv ]; then
  touch $HD/.done
  (flock 9; python3 - $ARM $SET $COV confidence $HD/happy.summary.csv >> $B/postproc/results/summary_happy.tsv <<'PY'
import csv, sys
arm, st, cov, sset, f = sys.argv[1:6]
cols = ["Type","Filter","TRUTH.TOTAL","TRUTH.TP","TRUTH.FN","QUERY.TOTAL","QUERY.FP","QUERY.UNK","FP.gt","FP.al","METRIC.Recall","METRIC.Precision","METRIC.Frac_NA","METRIC.F1_Score"]
for r in csv.DictReader(open(f)):
    if r["Filter"] == "PASS":
        print("\t".join([arm, st, cov, sset] + [r.get(c, "") for c in cols]))
PY
  ) 9>$B/postproc/results/summary_happy.lock
  rm -f $HD/harm.vcf.gz*; echo "== done $(date)"
else echo "hap.py FAILED, see $HD/happy.log"; exit 1; fi
