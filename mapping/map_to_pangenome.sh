#!/bin/bash

# ==============================================================================
# Mapping Pipeline Script
# ==============================================================================

# Resolve where this script lives so it works regardless of the working dir.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

# Load site-specific settings (reference locations, tool jars, how to load your
# environment). Copy config.example.sh -> config.sh in the repo root and edit it.
# Anything set there can still be overridden by the CLI flags below.
if [ -f "$REPO_ROOT/config.sh" ]; then
    source "$REPO_ROOT/config.sh"
fi

# --- Default Constants (config.sh values win; otherwise these fallbacks apply) ---
DEFAULT_ROOT="${PANGENOME_ROOT:-$REPO_ROOT/references}"
DEFAULT_GRAPH_NAME="${GRAPH_NAME:-hprc-v1.1-mc-grch38}"
DEFAULT_HAPL_NAME="${HAPL_NAME:-hprc-v1.1-mc-grch38.hapl}"
DEFAULT_PATHS_NAME="${PATHS_NAME:-hprc-v1.1-mc-grch38.paths_list.txt}"
DEFAULT_PREFIX="${SEQ_PREFIX:-GRCh38#0#}"
DEFAULT_ABRA_JAR="${ABRA2_JAR:-$DEFAULT_ROOT/abra2-2.23.jar}"
RENAME_SCRIPT="${RENAME_SCRIPT:-$SCRIPT_DIR/rename_bam_stream.py}"
DEFAULT_KMER=15
DEFAULT_WINDOW=5
SAMPLING_KMER=29
SAMPLING_WINDOW=11

# --- Usage Function ---
usage() {
    echo "Usage: $0 -S <sample> --runs <runs.tsv> -t <threads> [options]"
    echo "   or: $0 -s1 <fastq1> [-s2 <fastq2>] [-S <sample>] -t <threads> [options]"
    echo ""
    echo "Required:"
    echo "  -S   <string>  Sample name (required with --runs)"
    echo "  --runs <path>  TSV with header: sample, run, fq1, fq2 (fq2 empty for single-end)"
    echo "  -t   <int>     Number of threads"
    echo ""
    echo "Single run (instead of --runs):"
    echo "  -s1  <path>    Input Fastq R1"
    echo "  -s2  <path>    Input Fastq R2 (paired-end)"
    echo ""
    echo "  --keep-prerealn  Keep the merged BAM from before indel realignment"
    echo ""
    echo "Mapping minimizer (haplotype sampling is fixed at k=$SAMPLING_KMER, w=$SAMPLING_WINDOW):"
    echo "  -k   <int>     Minimizer k-mer length (default: $DEFAULT_KMER)"
    echo "  -w   <int>     Minimizer window length (default: $DEFAULT_WINDOW)"
    echo ""
    echo "Optional Inputs:"
    echo "  -r   <path>    Root directory for reference/graph files (default: $DEFAULT_ROOT)"
    echo "  -ref <path>    Reference FASTA file (default: looks for hg38.fa in root/ref)"
    echo ""
    echo "Optional Graph Files (IMPORTANT: These 3 files are a set. If you change -g, update -h, -p and -pre too):"
    echo "  -g   <name>    Graph GBZ filename (default: $DEFAULT_GRAPH_NAME.gbz)"
    echo "  -h   <name>    Haplotype filename (default: $DEFAULT_HAPL_NAME)"
    echo "  -p   <name>    Paths List filename (default: $DEFAULT_PATHS_NAME)"
    echo "  -pre <string>  Sequence prefix for renaming (default: $DEFAULT_PREFIX)"
    echo ""
    exit 1
}

kmer=$DEFAULT_KMER
window=$DEFAULT_WINDOW

# --- Parse command line arguments ---
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -s1)  s1="$2"; shift ;;
        -s2)  s2="$2"; shift ;;
        -S)   sample_base="$2"; shift ;;
        --runs) runs_tsv="$2"; shift ;;
        --keep-prerealn) keep_prerealn=1 ;;
        -tag)  tag="$2"; shift ;;
        -t)   threads="$2"; shift ;;
        -r)   r_in="$2"; shift ;;
        -g)   g_in="$2"; shift ;;
        -h)   h_in="$2"; shift ;;
        -p)   p_in="$2"; shift ;;
        -pre) pre_in="$2"; shift ;;
        -ref) ref_in="$2"; shift ;; 
        -k)   kmer="$2"; shift ;;   
        -w)   window="$2"; shift ;;  
        *) usage ;; 
    esac
    shift
done

# Check required inputs
if [ -z "$threads" ] || { [ -z "$s1" ] && [ -z "$runs_tsv" ]; }; then usage; fi
if [ -n "$runs_tsv" ] && [ -z "$sample_base" ]; then usage; fi

# --- Final Variable Assignment ---
root_dir=${r_in:-$DEFAULT_ROOT}
seqn_prefix=${pre_in:-$DEFAULT_PREFIX}
ABRA2_JAR=${DEFAULT_ABRA_JAR}

# --- INTELLIGENT PATH HANDLING ---

# 1. Handling Reference Genome (order: -ref flag > config REF_FA > files in root)
if [ -n "$ref_in" ]; then
    ref_fa="$ref_in"
elif [ -n "$REF_FA" ] && [ -f "$REF_FA" ]; then
    ref_fa="$REF_FA"
elif [ -f "$root_dir/hg38.fa" ]; then
    ref_fa="$root_dir/hg38.fa"
elif [ -f "$root_dir/genome.fa" ]; then
    # Auto-detect generic genome in root (common in test data)
    ref_fa="$root_dir/genome.fa"
else
    # Default Fallback
    ref_fa="$root_dir/ref/hg38.fa"
fi
ref_idx="${ref_fa}.fai"

# 2. Handling Graph Files (Fixes "File not found at .../ref/..." error)
gbz_name=${g_in:-"$DEFAULT_GRAPH_NAME.gbz"}
hapl_name=${h_in:-"$DEFAULT_HAPL_NAME"}
paths_name=${p_in:-"$DEFAULT_PATHS_NAME"}

# Function to check root first, then root/ref
find_file() {
    local fname=$1
    if [ -f "$root_dir/$fname" ]; then
        echo "$root_dir/$fname"
    elif [ -f "$root_dir/ref/$fname" ]; then
        echo "$root_dir/ref/$fname"
    else
        echo "MISSING"
    fi
}

gbz_path=$(find_file "$gbz_name")
hapl_path=$(find_file "$hapl_name")
ref_paths_list=$(find_file "$paths_name")

RUN_IDS=(); RUN_FQ1=(); RUN_FQ2=()
if [ -n "$runs_tsv" ]; then
    if [ ! -f "$runs_tsv" ]; then echo "ERROR: Runs TSV not found ($runs_tsv)"; exit 1; fi
    while IFS=$'\t' read -r t_sample t_run t_fq1 t_fq2; do
        [ "$t_sample" = "$sample_base" ] || continue
        RUN_IDS+=("$t_run"); RUN_FQ1+=("$t_fq1"); RUN_FQ2+=("$t_fq2")
    done < <(tail -n +2 "$runs_tsv" | tr -d '\r')
    if [ ${#RUN_IDS[@]} -eq 0 ]; then echo "ERROR: No runs for sample $sample_base in $runs_tsv"; exit 1; fi
else
    [ -n "$sample_base" ] || sample_base=$(basename "$s1" | sed -E 's/_R[12].*//; s/\.fastq.*//; s/\.fq.*//')
    if [ -z "$sample_base" ]; then
        echo "Error: Failed to extract sample name from $s1"; exit 1;
    fi
    RUN_IDS=("run1"); RUN_FQ1=("$s1"); RUN_FQ2=("$s2")
fi

sample="${sample_base}"
MAPTAG="k${kmer}_w${window}"
RES_DIR="results/${sample}"
RUNS_DIR="$RES_DIR/runs"

mkdir -p "$RES_DIR" "$RUNS_DIR"
mkdir -p "$RES_DIR/temp_kmc"

# Define Final Output Targets
KMC_OUT="$RES_DIR/${sample}.kff"
HAP_OUT="$RES_DIR/${sample}.haplotypes.gbz"
DIST_OUT="$RES_DIR/${sample}.dist"
MIN_OUT="$RES_DIR/${sample}.${MAPTAG}.min"
ZIPCODES="$RES_DIR/${sample}.${MAPTAG}.zipcodes"
SURJ_BAM="$RES_DIR/${sample}.${MAPTAG}.surj.bam"
FINAL_BAM="$RES_DIR/${sample}.${MAPTAG}.surj_realn.bam"

# --- Sanity Check Files ---
if [ "$gbz_path" == "MISSING" ]; then echo "ERROR: GBZ file not found ($gbz_name)"; exit 1; fi
if [ "$hapl_path" == "MISSING" ]; then echo "ERROR: Haplotype file not found ($hapl_name)"; exit 1; fi
if [ "$ref_paths_list" == "MISSING" ]; then echo "ERROR: Paths list not found ($paths_name)"; exit 1; fi
for f in "${RUN_FQ1[@]}" "${RUN_FQ2[@]}"; do
    if [ -n "$f" ] && [ ! -f "$f" ]; then echo "ERROR: Input Fastq not found ($f)"; exit 1; fi
done

# --- Environment Setup ---
# Define load_environment() in config.sh to make the required tools available
# (module loads, conda activate, PATH tweaks, ...). If it isn't defined we fall
# back to the LISC module setup this pipeline was developed on.
if declare -f load_environment > /dev/null; then
    load_environment
else
    ml load KMC vg Conda freebayes SAMtools BEDTools
    conda activate GATK-3.8.0
    conda activate freebayes
fi
echo "-------------------------------------------------------"
echo "Sample:         $sample"
echo "Reference:      $ref_fa"
echo "Graph:          $gbz_path"
echo "Output Dir:     $RES_DIR"
echo "-------------------------------------------------------"

# --- STEP 1: Count k-mers ---
if [ -f "$KMC_OUT" ] || [ -f "$HAP_OUT" ]; then
    echo "[SKIP] Step 1: KMC output already exists."
else
    echo "[RUN] Step 1: Counting k-mers"
    printf '%s\n' "${RUN_FQ1[@]}" "${RUN_FQ2[@]}" | grep -v '^$' > kmc_filelist.txt
    kmc -k${SAMPLING_KMER} -m2 -okff -t${threads} @kmc_filelist.txt "$RES_DIR/${sample}" "$RES_DIR/temp_kmc"
    rm -rf "$RES_DIR/temp_kmc" kmc_filelist.txt
fi

# --- STEP 2: Determine haplotypes ---
if [ -f "$HAP_OUT" ]; then
    echo "[SKIP] Step 2: Haplotype GBZ already exists."
else
    echo "[RUN] Step 2: Determining haplotypes"
    vg haplotypes -v 2 -t ${threads} \
        --num-haplotypes 32 --present-discount 0.9 --het-adjustment 0.05 \
        --absent-score 0.8 --include-reference --diploid-sampling --kmer-length ${SAMPLING_KMER} --window-leng ${SAMPLING_WINDOW} \
        -i ${hapl_path} -k ${KMC_OUT} -g ${HAP_OUT} ${gbz_path}
    [[ $? -eq 0 ]] && rm -f "${KMC_OUT}"
fi

# --- STEP 3: Index graph (distance) ---
if [ -f "$DIST_OUT" ]; then
    echo "[SKIP] Step 3: Distance index exists."
else
    vg index -j ${DIST_OUT} ${HAP_OUT}
fi

# --- STEP 4: Index graph (minimizer) ---
if [ -f "$MIN_OUT" ] || [ -f "$SURJ_BAM" ] || [ -f "$FINAL_BAM" ]; then
    echo "[SKIP] Step 4: Minimizer index exists."
else
    echo "[RUN] Step 4: Creating Minimizer Index"
    vg minimizer -p -t ${threads} -d ${DIST_OUT} -o ${MIN_OUT} -z ${ZIPCODES} ${HAP_OUT} -k ${kmer} -w ${window}
fi

# --- STEP 5: Mapping with vg giraffe (one GAF per run) ---
GIRAFFE_OPTS="--progress --sample ${sample} --output-format gaf -Z ${HAP_OUT} -d ${DIST_OUT} -m ${MIN_OUT} -z ${ZIPCODES} -t ${threads}"
for i in "${!RUN_IDS[@]}"; do
    run="${RUN_IDS[$i]}"
    GAF_OUT="$RUNS_DIR/${run}.${MAPTAG}.gaf.gz"
    RUN_BAM="$RUNS_DIR/${run}.${MAPTAG}.surj.bam"
    if [ -f "$GAF_OUT" ] || [ -f "$RUN_BAM" ] || [ -f "$SURJ_BAM" ] || [ -f "$FINAL_BAM" ]; then
        echo "[SKIP] Step 5: GAF for run $run exists."
    else
        echo "[RUN] Step 5: Mapping run $run with vg giraffe"
        if [ -n "${RUN_FQ2[$i]}" ]; then
            vg giraffe $GIRAFFE_OPTS -f "${RUN_FQ1[$i]}" -f "${RUN_FQ2[$i]}" | gzip > "${GAF_OUT}"
        else
            vg giraffe $GIRAFFE_OPTS -f "${RUN_FQ1[$i]}" | gzip > "${GAF_OUT}"
        fi
    fi
done

# --- STEP 6: Surject (one GAF to one BAM per run), then merge per sample ---
if [ -f "$SURJ_BAM" ] || [ -f "$FINAL_BAM" ]; then
    echo "[SKIP] Step 6: Surjected BAM exists."
else
    # SAFETY CHECK: Ensure bamleftalign is actually runnable
    if ! command -v bamleftalign &> /dev/null; then
        echo "CRITICAL ERROR: 'bamleftalign' command not found."
        echo "Please ensure freebayes is installed in your current Conda environment."
        exit 1
    fi

    RUN_BAMS=()
    for i in "${!RUN_IDS[@]}"; do
        run="${RUN_IDS[$i]}"
        GAF_OUT="$RUNS_DIR/${run}.${MAPTAG}.gaf.gz"
        RUN_BAM="$RUNS_DIR/${run}.${MAPTAG}.surj.bam"
        RUN_BAMS+=("$RUN_BAM")
        if [ -f "$RUN_BAM" ]; then
            echo "[SKIP] Step 6: BAM for run $run exists."
            continue
        fi

        vg surject \
            -F "${ref_paths_list}" -x ${HAP_OUT} -t ${threads} \
            --sam-output --gaf-input --sample ${sample} \
            --read-group "ID:${run} LB:${run} SM:${sample} PL:illumina" \
            --prune-low-cplx ${GAF_OUT} > "$RES_DIR/tmp.bam"

        python3 "$RENAME_SCRIPT" -f "${ref_idx}" -p "$seqn_prefix" < "$RES_DIR/tmp.bam" > "$RES_DIR/renamed.bam"

        bamleftalign --fasta-reference "${ref_fa}" --compressed < "$RES_DIR/renamed.bam" | \
        samtools sort -m 1G --threads ${threads} -o ${RUN_BAM}

        if [ $? -ne 0 ]; then
            echo "Error: bamleftalign or samtools sort failed."
            exit 1
        fi
    done

    samtools merge -f --threads ${threads} ${SURJ_BAM} "${RUN_BAMS[@]}"
    samtools index ${SURJ_BAM}
    rm -f "$RES_DIR/tmp.bam" "$RES_DIR/renamed.bam" "${RUN_BAMS[@]}" ${MIN_OUT} ${ZIPCODES}
fi

# --- STEP 8 & 9: GATK + ABRA2 Realignment ---
TARGETS_BED="$RES_DIR/${sample}.realn_targets.bed"

if [ -f "$FINAL_BAM" ]; then
    echo "[SKIP] Steps 8-9: Final BAM exists."
else
    echo "[RUN] Step 8: GATK RealignerTargetCreator"
    
    # GATK 3.8 (RealignerTargetCreator) needs Java 8. Set JAVA8 and GATK3_JAR in
    # config.sh; the fallbacks below point at the LISC conda env used originally.
    JAVA8_CMD="${JAVA8:-/lisc/opt/sw/software/Conda/Miniforge3/envs/GATK-3.8.0/bin/java}"
    GATK_JAR="${GATK3_JAR:-/lisc/opt/sw/software/Conda/Miniforge3/envs/GATK-3.8.0/opt/gatk-3.8/GenomeAnalysisTK.jar}"
    if [ ! -f "$GATK_JAR" ]; then echo "ERROR: GATK Jar not found at $GATK_JAR"; exit 1; fi

    $JAVA8_CMD -Xmx16G  -Djava.io.tmpdir="$RES_DIR" -jar "$GATK_JAR" -T RealignerTargetCreator \
        -nt ${threads} -R "${ref_fa}" -I ${SURJ_BAM} --out "$RES_DIR/intervals"
    if [ ! -f "$RES_DIR/intervals" ]; then
        echo "ERROR: GATK failed to create intervals file. Check logs."
        exit 1
    fi

    awk -F '[:-]' 'BEGIN { OFS = "\t" } { if( $3 == "") { print $1, $2-1, $2 } else { print $1, $2-1, $3}}' "$RES_DIR/intervals" > "$RES_DIR/tmp.bed"
    bedtools slop -i "$RES_DIR/tmp.bed" -g "${ref_idx}" -b "160" > ${TARGETS_BED}

    echo "[RUN] Step 9: ABRA2"
    # Handling Java version switch if necessary (ABRA2 works on Java 8 usually)
    # ml unload java && ml load java # Uncomment if your system needs this
    
    java -Xmx26G -jar "${ABRA2_JAR}" --targets ${TARGETS_BED} --in ${SURJ_BAM} --out ${FINAL_BAM} --ref "${ref_fa}" --threads ${threads}
    
    samtools index ${FINAL_BAM}
    [ -n "$keep_prerealn" ] || rm -f ${SURJ_BAM} "${SURJ_BAM}.bai"
    rm -f ${TARGETS_BED} "$RES_DIR/intervals" "$RES_DIR/tmp.bed"
fi
echo "[RUN] Step Stats"
# --- Step 1: Generate Detailed Stats ---
STATS_OUT="$RES_DIR/${sample}_samtools.stats"
FLAG_OUT="$RES_DIR/${sample}_flagstat.txt"

samtools stats -@ "$threads" "$FINAL_BAM" > "$STATS_OUT"
# flagstat is good for a quick sanity check but redundant if using stats
samtools flagstat "$FINAL_BAM" > "$FLAG_OUT"

# --- Step 2: Extract Summary Metrics ---
# Using 'awk' here is safer than 'cut' because it handles variable whitespace
total_reads=$(grep "^SN	raw total sequences:" "$STATS_OUT" | awk -F'\t' '{print $3}')
mapped_reads=$(grep "^SN	reads mapped:" "$STATS_OUT" | awk -F'\t' '{print $3}')
unmapped_reads=$(grep "^SN	reads unmapped:" "$STATS_OUT" | awk -F'\t' '{print $3}')

# --- Step 3: Extract High Quality (MAPQ >= 30) ---
# We use 'samtools view -c' (count) with '-q 30' (min quality)
# This is the most reliable way to get this number.
gt30_reads=$(samtools view -@ "$threads" -c -q 30 "$FINAL_BAM")

# --- Step 4: Extract MAPQ Distribution for Plotting ---
# This creates your 2-column TSV for the heatmap/line plots later
grep "^MQ" "$STATS_OUT" | awk -F'\t' '{print $2"\t"$3}' > "$RES_DIR/${sample}_mapq_distribution.tsv"

# --- Step 5: Save the Summary Row ---
# Formatting with column headers (optional, but helpful for debugging)
echo -e "Kmer\tWindow\tTotal\tMapped\tUnmapped\tMQ30" > "$RES_DIR/row_k${kmer}_w${window}.tsv"
echo -e "${kmer}\t${window}\t${total_reads}\t${mapped_reads}\t${unmapped_reads}\t${gt30_reads}" >> "$RES_DIR/row_k${kmer}_w${window}.tsv"

# --- Step 6: The "Nuclear" Cleanup ---
# CRITICAL: We do this AFTER we've counted gt30_reads
echo "[RUN] Cleaning up intermediate files in $RES_DIR..."

# Delete every file that IS NOT a .txt or .tsv (Removes BAM, BAI, GAF, etc.)
#find "$RES_DIR" -type f ! -name "*.txt" ! -name "*.tsv" -delete

# Remove empty subdirectories if they exist
#find "$RES_DIR" -mindepth 1 -type d -empty -delete
