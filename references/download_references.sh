#!/bin/bash
# ==============================================================================
# download_references.sh
# ------------------------------------------------------------------------------
# Downloads the large, PUBLIC reference files that map_to_pangenome.sh needs.
# These files are NOT stored in git (they are tens of GB). Run this once after
# cloning the repo. Everything lands in this "references/" directory.
#
#   Downloads (~30 GB total):
#     * hg38.fa                      linear reference genome (UCSC)   ~3 GB
#     * hg38.fa.fai / hg38.dict      indexes for the above (built here)
#     * hprc-v1.1-mc-grch38.gbz      HPRC pangenome graph             ~3 GB
#     * hprc-v1.1-mc-grch38.hapl     haplotype-sampling index        ~20 GB
#
#   (hprc-v1.1-mc-grch38.paths_list.txt is small and ships with the repo.)
#
# Requires: wget (or curl), gzip, samtools. samtools is used to build the .fai
# and .dict; if it is not on your PATH, load it first (e.g. `ml load SAMtools`).
# ==============================================================================
set -euo pipefail

# Download into the directory this script lives in.
DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DEST"

HPRC_BASE="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.1-mc-grch38"
HG38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"

# Prefer wget, fall back to curl.
fetch() {  # fetch <url> <output>
    local url="$1" out="$2"
    if [ -f "$out" ]; then echo "[SKIP] $out already exists"; return; fi
    echo "[GET ] $out"
    if command -v wget >/dev/null 2>&1; then
        wget -c -O "$out" "$url"
    else
        curl -L -C - -o "$out" "$url"
    fi
}

echo "==> Reference genome (hg38)"
fetch "$HG38_URL" "hg38.fa.gz"
if [ ! -f "hg38.fa" ]; then
    echo "[GUNZIP] hg38.fa.gz"
    gunzip -k hg38.fa.gz
fi

echo "==> Indexing hg38 (.fai + .dict)"
command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not on PATH (try: ml load SAMtools)"; exit 1; }
[ -f "hg38.fa.fai" ] || samtools faidx hg38.fa
[ -f "hg38.dict" ]   || samtools dict hg38.fa -o hg38.dict

echo "==> HPRC pangenome graph"
fetch "$HPRC_BASE/hprc-v1.1-mc-grch38.gbz"  "hprc-v1.1-mc-grch38.gbz"
fetch "$HPRC_BASE/hprc-v1.1-mc-grch38.hapl" "hprc-v1.1-mc-grch38.hapl"

echo ""
echo "Done. Reference files are in: $DEST"
echo "If you use a config.sh, point PANGENOME_ROOT at this directory."
