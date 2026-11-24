#!/bin/bash

# ==============================================================================
# Mapping pipeline: reads -> personalized HPRC graph -> BAM on GRCh38 coordinates
#
#   1. kmc            count k-mers over all runs of the sample (k=29)
#   2. vg haplotypes  sample one personalized graph per sample
#   3. vg index       distance index
#   4. vg minimizer   minimizer index (mapping k/w)
#   5. vg giraffe     map each run, write BAM directly (-o BAM, surjected on the
#                     fly), strip the contig prefix, sort; merge runs per sample
#   6. samtools stats summary
# ==============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# Site-specific settings live in config.sh (git-ignored). See config.example.sh.
if [ -f "$REPO_ROOT/config.sh" ]; then
    source "$REPO_ROOT/config.sh"
fi

# --- Defaults (config.sh values win; CLI flags win over config.sh) ---
DEFAULT_ROOT="${PANGENOME_ROOT:-$REPO_ROOT/references}"
DEFAULT_GRAPH_NAME="${GRAPH_NAME:-hprc-v1.1-mc-grch38}"
DEFAULT_HAPL_NAME="${HAPL_NAME:-hprc-v1.1-mc-grch38.hapl}"
DEFAULT_PATHS_NAME="${PATHS_NAME:-hprc-v1.1-mc-grch38.paths_list.txt}"
DEFAULT_PREFIX="${SEQ_PREFIX:-GRCh38#0#}"
RENAME_SCRIPT="${RENAME_SCRIPT:-$SCRIPT_DIR/rename_bam_stream.py}"
# Mapping minimizer: vg giraffe's own defaults unless set in config.sh / -k -w.
DEFAULT_KMER="${MAP_KMER:-29}"
DEFAULT_WINDOW="${MAP_WINDOW:-11}"
# Haplotype sampling (fixed k/w; number of candidate haplotypes configurable).
SAMPLING_KMER=29
SAMPLING_WINDOW=11
DEFAULT_NUM_HAPS="${NUM_HAPLOTYPES:-4}"

usage() {
    echo "Usage: $0 -S <sample> --runs <runs.tsv> -t <threads> [options]"
    echo ""
    echo "Required:"
    echo "  -S   <string>  Sample name"
    echo "  --runs <path>  TSV with header: sample, run, fq1, fq2 (fq2 empty for single-end)"
    echo "  -t   <int>     Number of threads"
    echo ""
    echo "Mapping / sampling:"
    echo "  -k   <int>     Mapping minimizer k (default: $DEFAULT_KMER; config MAP_KMER)"
    echo "  -w   <int>     Mapping minimizer w (default: $DEFAULT_WINDOW; config MAP_WINDOW)"
    echo "  --num-haplotypes <int>  Candidate haplotypes for sampling (default: $DEFAULT_NUM_HAPS)"
    echo "                 (haplotype sampling itself is fixed at k=$SAMPLING_KMER, w=$SAMPLING_WINDOW)"
    echo ""
    echo "Optional inputs:"
    echo "  -r   <path>    Root directory for graph files (default: $DEFAULT_ROOT)"
    echo "  -ref <path>    Reference FASTA (default: \$REF_FA, else root/hg38.fa)"
    echo "  -g   <name>    Graph GBZ filename (default: $DEFAULT_GRAPH_NAME.gbz)"
    echo "  -h   <name>    Haplotype index filename (default: $DEFAULT_HAPL_NAME)"
    echo "  -p   <name>    Paths list filename (default: $DEFAULT_PATHS_NAME)"
    echo "  -pre <string>  Contig prefix to strip (default: $DEFAULT_PREFIX)"
    echo "  (-g, -h, -p and -pre are a set: change them together)"
    exit 1
}

kmer=$DEFAULT_KMER
window=$DEFAULT_WINDOW
num_haps=$DEFAULT_NUM_HAPS
sample=""; runs_tsv=""; threads=""
r_in=""; g_in=""; h_in=""; p_in=""; pre_in=""; ref_in=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -S)     sample="$2"; shift ;;
        --runs) runs_tsv="$2"; shift ;;
        -t)     threads="$2"; shift ;;
        -k)     kmer="$2"; shift ;;
        -w)     window="$2"; shift ;;
        --num-haplotypes) num_haps="$2"; shift ;;
        -r)     r_in="$2"; shift ;;
        -g)     g_in="$2"; shift ;;
        -h)     h_in="$2"; shift ;;
        -p)     p_in="$2"; shift ;;
        -pre)   pre_in="$2"; shift ;;
        -ref)   ref_in="$2"; shift ;;
        *) usage ;;
    esac
    shift
done
if [ -z "$threads" ] || [ -z "$runs_tsv" ] || [ -z "$sample" ]; then usage; fi

root_dir=${r_in:-$DEFAULT_ROOT}
seqn_prefix=${pre_in:-$DEFAULT_PREFIX}

# --- Reference FASTA: -ref > config REF_FA > files in root ---
if [ -n "$ref_in" ]; then ref_fa="$ref_in"
elif [ -n "${REF_FA:-}" ] && [ -f "$REF_FA" ]; then ref_fa="$REF_FA"
elif [ -f "$root_dir/hg38.fa" ]; then ref_fa="$root_dir/hg38.fa"
elif [ -f "$root_dir/genome.fa" ]; then ref_fa="$root_dir/genome.fa"
else ref_fa="$root_dir/ref/hg38.fa"
fi
ref_idx="${ref_fa}.fai"

# --- Graph files: look in root, then root/ref ---
find_file() {
    if [ -f "$root_dir/$1" ]; then echo "$root_dir/$1"
    elif [ -f "$root_dir/ref/$1" ]; then echo "$root_dir/ref/$1"
    else echo "MISSING"; fi
}
gbz_path=$(find_file "${g_in:-$DEFAULT_GRAPH_NAME.gbz}")
hapl_path=$(find_file "${h_in:-$DEFAULT_HAPL_NAME}")
ref_paths_list=$(find_file "${p_in:-$DEFAULT_PATHS_NAME}")
for f in "$gbz_path" "$hapl_path" "$ref_paths_list"; do
    [ "$f" = "MISSING" ] && { echo "ERROR: graph/haplotype/paths file not found under $root_dir"; exit 1; }
done
for f in "$ref_fa" "$ref_idx" "$runs_tsv"; do
    [ -f "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# --- Read this sample's runs from the TSV ---
RUN_IDS=(); RUN_FQ1=(); RUN_FQ2=()
while IFS=$'\t' read -r t_sample t_run t_fq1 t_fq2; do
    [ "$t_sample" = "$sample" ] || continue
    RUN_IDS+=("$t_run"); RUN_FQ1+=("$t_fq1"); RUN_FQ2+=("$t_fq2")
done < <(tail -n +2 "$runs_tsv" | tr -d '\r')
[ ${#RUN_IDS[@]} -gt 0 ] || { echo "ERROR: no runs for $sample in $runs_tsv"; exit 1; }
for f in "${RUN_FQ1[@]}" "${RUN_FQ2[@]}"; do
    [ -z "$f" ] || [ -f "$f" ] || { echo "ERROR: fastq not found: $f"; exit 1; }
done

MAPTAG="k${kmer}_w${window}"
RES_DIR="results/${sample}"
RUNS_DIR="$RES_DIR/runs"
mkdir -p "$RES_DIR" "$RUNS_DIR" "$RES_DIR/temp_kmc"

KMC_OUT="$RES_DIR/${sample}.kff"
HAP_OUT="$RES_DIR/${sample}.haplotypes.gbz"
DIST_OUT="$RES_DIR/${sample}.dist"
MIN_OUT="$RES_DIR/${sample}.${MAPTAG}.min"
ZIPCODES="$RES_DIR/${sample}.${MAPTAG}.zipcodes"
FINAL_BAM="$RES_DIR/${sample}.${MAPTAG}.bam"

# --- Environment ---
if declare -f load_environment > /dev/null; then
    load_environment
else
    ml load KMC vg SAMtools
fi
echo "Sample: $sample | graph: $gbz_path | haplotypes: $num_haps | mapping k/w: $kmer/$window"
date

# --- Step 1: k-mer counts over all runs ---
if [ -f "$KMC_OUT" ] || [ -f "$HAP_OUT" ]; then
    echo "[SKIP] Step 1: kmc"
else
    echo "[RUN] Step 1: kmc"
    printf '%s\n' "${RUN_FQ1[@]}" "${RUN_FQ2[@]}" | grep -v '^$' > "$RES_DIR/kmc_filelist.txt"
    kmc -k${SAMPLING_KMER} -m2 -okff -t${threads} @"$RES_DIR/kmc_filelist.txt" "$RES_DIR/${sample}" "$RES_DIR/temp_kmc" || exit 1
    rm -rf "$RES_DIR/temp_kmc" "$RES_DIR/kmc_filelist.txt"
fi

# --- Step 2: personalized graph ---
if [ -f "$HAP_OUT" ]; then
    echo "[SKIP] Step 2: haplotypes"
else
    echo "[RUN] Step 2: vg haplotypes (--num-haplotypes $num_haps)"
    vg haplotypes -v 2 -t ${threads} \
        --num-haplotypes ${num_haps} --present-discount 0.9 --het-adjustment 0.05 \
        --absent-score 0.8 --include-reference --diploid-sampling \
        --kmer-length ${SAMPLING_KMER} --window-length ${SAMPLING_WINDOW} \
        -i ${hapl_path} -k ${KMC_OUT} -g ${HAP_OUT} ${gbz_path} || exit 1
    rm -f "${KMC_OUT}"
fi

# --- Step 3: distance index ---
if [ -f "$DIST_OUT" ]; then
    echo "[SKIP] Step 3: distance index"
else
    echo "[RUN] Step 3: distance index"
    vg index -j ${DIST_OUT} ${HAP_OUT} || exit 1
fi

# --- Step 4: minimizer index ---
if [ -f "$MIN_OUT" ] || [ -f "$FINAL_BAM" ]; then
    echo "[SKIP] Step 4: minimizer"
else
    echo "[RUN] Step 4: minimizer k=$kmer w=$window"
    vg minimizer -p -t ${threads} -d ${DIST_OUT} -o ${MIN_OUT} -z ${ZIPCODES} ${HAP_OUT} -k ${kmer} -w ${window} || exit 1
fi

# --- Step 5: giraffe -> BAM, one per run, then merge ---
if [ -f "$FINAL_BAM" ]; then
    echo "[SKIP] Step 5: BAM exists"
else
    RUN_BAMS=()
    for i in "${!RUN_IDS[@]}"; do
        run="${RUN_IDS[$i]}"
        RUN_BAM="$RUNS_DIR/${run}.${MAPTAG}.bam"
        RUN_BAMS+=("$RUN_BAM")
        [ -f "$RUN_BAM" ] && { echo "[SKIP] run $run"; continue; }
        echo "[RUN] Step 5: vg giraffe -o BAM, run $run"
        FQARGS=(-f "${RUN_FQ1[$i]}")
        [ -n "${RUN_FQ2[$i]}" ] && FQARGS+=(-f "${RUN_FQ2[$i]}")
        vg giraffe --progress --sample ${sample} --read-group ${run} -o BAM \
            -Z ${HAP_OUT} -d ${DIST_OUT} -m ${MIN_OUT} -z ${ZIPCODES} \
            --ref-paths ${ref_paths_list} --prune-low-cplx -t ${threads} "${FQARGS[@]}" \
          | samtools view -h - \
          | python3 "$RENAME_SCRIPT" -f "${ref_idx}" -p "$seqn_prefix" \
          | samtools sort -m 1G --threads ${threads} -o "${RUN_BAM}.tmp.bam" -
        rc=("${PIPESTATUS[@]}")
        if [ "${rc[0]}" -ne 0 ] || [ "${rc[2]}" -ne 0 ] || [ "${rc[3]}" -ne 0 ] || [ ! -s "${RUN_BAM}.tmp.bam" ]; then
            echo "ERROR: giraffe/rename/sort failed for run $run (exit codes: ${rc[*]})"
            rm -f "${RUN_BAM}.tmp.bam"; exit 1
        fi
        mv "${RUN_BAM}.tmp.bam" "${RUN_BAM}"
    done
    if [ ${#RUN_BAMS[@]} -eq 1 ]; then
        mv "${RUN_BAMS[0]}" "${FINAL_BAM}"
    else
        samtools merge -f --threads ${threads} ${FINAL_BAM} "${RUN_BAMS[@]}" && rm -f "${RUN_BAMS[@]}" || exit 1
    fi
    samtools index ${FINAL_BAM}
    rm -f ${MIN_OUT} ${ZIPCODES}
fi

# --- Step 6: summary statistics ---
echo "[RUN] Step 6: stats"
STATS_OUT="$RES_DIR/${sample}_samtools.stats"
samtools stats -@ "$threads" "$FINAL_BAM" > "$STATS_OUT"
samtools flagstat "$FINAL_BAM" > "$RES_DIR/${sample}_flagstat.txt"
total_reads=$(grep "^SN	raw total sequences:" "$STATS_OUT" | awk -F'\t' '{print $3}')
mapped_reads=$(grep "^SN	reads mapped:" "$STATS_OUT" | awk -F'\t' '{print $3}')
unmapped_reads=$(grep "^SN	reads unmapped:" "$STATS_OUT" | awk -F'\t' '{print $3}')
gt30_reads=$(samtools view -@ "$threads" -c -q 30 "$FINAL_BAM")
grep "^MQ" "$STATS_OUT" | awk -F'\t' '{print $2"\t"$3}' > "$RES_DIR/${sample}_mapq_distribution.tsv"
echo -e "Kmer\tWindow\tTotal\tMapped\tUnmapped\tMQ30" > "$RES_DIR/row_${MAPTAG}.tsv"
echo -e "${kmer}\t${window}\t${total_reads}\t${mapped_reads}\t${unmapped_reads}\t${gt30_reads}" >> "$RES_DIR/row_${MAPTAG}.tsv"
echo "done"; date
