#!/bin/bash
# Score one finished set with hap.py (vcfeval engine, genotype-level, GIAB truth).
#
#   1. build the "confidence" call set from raw.vcf.gz:
#        FILTER PASS (or '.'), FORMAT/DP >= MIN_DP, het calls with allele balance 0.3-0.7
#   2. drop ALT=* records, run hap.py (rtg heap capped by RTGMEM)
#   3. append the PASS rows to $B/postproc/results/summary_happy.tsv
#
# Array: task N reads line N+1 of $LIST, each line = "<arm> <set>", set name ends in cov<COV>.
# MIN_DP is the minimum depth for a call to be kept (default 5). Change it per submission, e.g.
#   sbatch --array=0-N --export=ALL,LIST=list.txt,MIN_DP=1 happy_array.sbatch
# Results go to <set>/happy/confidence_DP<MIN_DP>/, so different thresholds never overwrite each other.
#SBATCH --cpus-per-task=8
#SBATCH --mem=100G
#SBATCH --time=08:00:00
set -uo pipefail
# Paths come from benchmarking/config.sh (copy config.example.sh); submit from benchmarking/ or set BENCH_CONFIG.
source "${BENCH_CONFIG:-${SLURM_SUBMIT_DIR:?submit from benchmarking/}/config.sh}"
MIN_DP=${MIN_DP:-5}
TAG=confidence_DP${MIN_DP}
read ARM SET < <(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" ${LIST:?LIST not set})
COV=${SET##*cov}; W=$B/postproc/results/$ARM/$SET; HD=$W/happy/$TAG
[ -e $W/.done_sets ] && [ -s $W/raw.vcf.gz ] || { echo "not finished: $ARM $SET"; exit 1; }
[ -e $HD/.done ] && { echo "already scored $ARM $SET ($TAG)"; exit 0; }
ml load BCFtools SAMtools >/dev/null 2>&1
mkdir -p $HD; rm -f $HD/happy.*; echo "== hap.py $ARM $SET DP>=$MIN_DP $(date)"
bcftools view -f PASS,. -i "FORMAT/DP>=$MIN_DP" $W/raw.vcf.gz -Ov | awk -f $BP/hetfilt_0.3_0.7.awk \
  | bcftools view -Oz -o $W/$TAG.vcf.gz && bcftools index -tf $W/$TAG.vcf.gz || exit 1
echo "records: $(bcftools view -H $W/$TAG.vcf.gz | wc -l)"
# drop GATK's ALT=* (spanning-deletion) records, which the truth set does not contain
bcftools view -e 'ALT="*"' $W/$TAG.vcf.gz -Oz -o $HD/harm.vcf.gz && bcftools index -tf $HD/harm.vcf.gz || exit 1
RTGMEM=64g apptainer exec --bind "$BIND_PATHS" --env RTG_TALKBACK=false $HAPPY_SIF \
  hap.py --engine vcfeval --engine-vcfeval-path $BP/rtg-wrapped-mem --engine-vcfeval-template $T_/hg38.sdf \
    --gender male -r $REF_FA -f $T_/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed \
    --stratification $T_/stratifications/GRCh38@all/GRCh38-all-stratifications.abspath.tsv \
    -o $HD/happy --threads $SLURM_CPUS_PER_TASK $T_/GRCh38_HG2-T2TQ100-V1.1_smvar.vcf.gz $HD/harm.vcf.gz > $HD/happy.log 2>&1
if [ -s $HD/happy.summary.csv ]; then
  touch $HD/.done
  (flock 9; python3 - $ARM $SET $COV $TAG $HD/happy.summary.csv >> $B/postproc/results/summary_happy.tsv <<'PY'
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
