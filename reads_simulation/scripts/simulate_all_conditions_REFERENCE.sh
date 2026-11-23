#!/bin/bash
# REFERENCE / DOCUMENTATION SCRIPT - the whole read simulation for all 9 conditions in one place, with the exact commands and parameters
# that were run (they were run as separate SLURM steps: 01a_fragments, 01_simulate, 02_trim_merge, 03_subsample; this script is not
# meant to be re-run as is: one condition takes many hours and ~100 GB of temporary files).
#
#   conditions = 3 read-length distributions (short, mid, long) x 3 damage levels (low, mid, high) = 9 simulations at 10x,
#   then nested sub-sampling to 7x, 5x, 1x.  Merged reads only are used for mapping.
set -euo pipefail

# ---------- inputs ----------
GENOME=/path/to/hg002v1.1.fasta
#   HG002 T2T-Q100 v1.1 assembly, 47 contigs, 6.0 Gb: maternal + paternal haplotypes + chrM, passed to gargammel as ONE haploid "endo" file
#   (copied to node-local disk first: reading it from the network disk was ~200x slower).
GARG=/path/to/gargammel
GENOME_BP=2999712359        # haploid genome length used for the coverage arithmetic
COV_MAX=10                  # simulate once at 10x, sub-sample for lower coverage
MINSIZE=30; MAXSIZE=291     # fragment length range (bp)
READLEN=125; PLATFORM=HS25  # paired-end 2 x 125 bp, ART Illumina HiSeq 2500 125 bp error profile
SEED_BASE=1000
ADAPTF=AGATCGGAAGAGCACACGTCTGAACTCCAGTCACCGATTCGATCTCGTATGCCGTCTTCTGCTTG
ADAPTR=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTT

# ---------- the 3 read-length conditions (fragSim length arguments) ----------
#   short: empirical distribution Ind1 (median 46 / mean 55 bp)   mid: empirical Ind2 (median 67 / mean 84 bp)
#   long : lognormal, meanlog 4.6 sdlog 0.4 (median ~100 / mean 107 bp after truncation to 30-291 bp)
declare -A LENARGS=( [short]="-f length/short.freq" [mid]="-f length/mid.freq" [long]="--loc 4.6 --scale 0.4" )
declare -A MEANLEN=( [short]=54.98 [mid]=84.22 [long]=107.06 )
# ---------- the 3 damage levels: Briggs model "-damage v,l,d,s" (double-stranded library), calibrated to 7 / 15 / 25 % C->T at the first base ----------
declare -A DAMAGE=( [low]="0.024,0.36,0.0031,0.206" [mid]="0.024,0.36,0.0066,0.447" [high]="0.024,0.36,0.0108,0.76" )

for COND in short mid long; do
  # number of fragments for 10x of the haploid genome length
  N=$(python3 -c "print(round($COV_MAX*$GENOME_BP/${MEANLEN[$COND]}))")      # short 545.6 M, mid 356.2 M, long 280.2 M fragments

  # 1) fragments: drawn from the assembly in proportion to contig length; ONE fragment set per length condition, shared by the 3 damage levels
  $GARG/src/fragSim --seed $((SEED_BASE+1)) -tag e -n $N -m $MINSIZE -M $MAXSIZE ${LENARGS[$COND]} $GENOME | gzip > $COND.frags.fa.gz

  for DAM in low mid high; do
    # 2) post-mortem damage on the fragments (deamination at the ends; C->T at 5', G->A at 3')
    $GARG/src/deamSim --seed $((SEED_BASE+2)) -damage ${DAMAGE[$DAM]} $COND.frags.fa.gz | gzip > $COND.$DAM.damaged.fa.gz

    # 3) add sequencing adapters so that fragments shorter than the read length read through into adapter
    $GARG/src/adptSim --seed $((SEED_BASE+3)) -f $ADAPTF -s $ADAPTR -l $READLEN -artp $COND.$DAM.adapt.fa.gz $COND.$DAM.damaged.fa.gz

    # 4) Illumina read simulation (default ART HS25 substitution AND small indel errors; quality scores not used)
    #    (conditions short/low, short/mid, mid/high, long/low, long/mid: this step was run in 20 parts with different seeds, 04_art_split/part/join.sh,
    #     because single ART runs stalled; parameters identical)
    $GARG/art_src_MountRainier/art_illumina -rs $((SEED_BASE+5)) -ss $PLATFORM -amp -na -p -i $COND.$DAM.adapt.fa.gz -l $READLEN -c 1 -qs 0 -qs2 0 -gzl 4 -o $COND.$DAM.art
    #    -> $COND.$DAM.art1.fq.gz / art2.fq.gz  (paired reads, 2 x 125 bp)

    # 5) adapter trimming and read merging (ancient-DNA mode): fragments up to ~239 bp can be merged (2 x 125 bp minus minimum overlap);
    #    merged reads  = <prefix>.fq.gz   (these are the reads used for mapping)   unmerged pairs = <prefix>_r1/_r2.fq.gz (not used)
    leeHom --ancientdna -fq1 $COND.$DAM.art1.fq.gz -fq2 $COND.$DAM.art2.fq.gz -fqo $COND.$DAM.leehom -t 4 --seed $((SEED_BASE)) --log $COND.$DAM.leehom.log

    # 6) coverage levels as nested random sub-samples of the 10x merged reads (same seed -> 1x inside 5x inside 7x inside 10x)
    mkdir -p final/$COND/$DAM/cov10; cp $COND.$DAM.leehom.fq.gz final/$COND/$DAM/cov10/merged.fq.gz
    for COV in 7 5 1; do
      mkdir -p final/$COND/$DAM/cov$COV
      seqtk sample -s 321 final/$COND/$DAM/cov10/merged.fq.gz $(python3 -c "print($COV/$COV_MAX)") | gzip > final/$COND/$DAM/cov$COV/merged.fq.gz
    done
  done
done
# Result: reads/final/<short|mid|long>/<low|mid|high>/cov<10|7|5|1>/merged.fq.gz  (36 single-end merged read sets used for the benchmark)
