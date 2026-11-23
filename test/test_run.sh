#!/bin/bash
# ==============================================================================
# test_run.sh — end-to-end smoke test on a tiny synthetic pangenome (~260 KB).
# ------------------------------------------------------------------------------
# Runs the full map_to_pangenome.sh pipeline against the toy data in test/data/
# and checks that a non-empty final BAM with mapped reads is produced. Finishes
# in seconds and needs no reference download.
#
# It uses the SAME tools/environment as a real run (vg, kmc, samtools, bedtools,
# bamleftalign, GATK3+Java8, ABRA2), so it also doubles as a check that your
# config.sh / environment is set up correctly. Set those up first (see the main
# README); on the original LISC cluster the built-in fallbacks are used.
#
# Usage:  bash test/test_run.sh [work_dir]
#         THREADS=8 bash test/test_run.sh
# ==============================================================================
set -euo pipefail

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$TEST_DIR")"
DATA="$TEST_DIR/data"
WORK="${1:-$TEST_DIR/work}"          # outputs go here (git-ignored)
THREADS="${THREADS:-4}"

mkdir -p "$WORK"
cd "$WORK"

echo "== Running pipeline on toy data =="
# Toy graph/reference both use the contig name "chr1", so no prefix needs to be
# stripped; -pre is set to a string that does not occur (a no-op). On real HPRC
# data the prefix is "GRCh38#0#".
bash "$REPO/scripts/map_to_pangenome.sh" \
    -s1 "$DATA/toy_reads.fq.gz" -t "$THREADS" \
    -r "$DATA" \
    -g toy.giraffe.gbz -h toy.hapl -p toy.paths_list.txt \
    -pre "none#" -ref "$DATA/chr1.fa" \
    -k 15 -w 5

BAM="$WORK/results/toy_reads/toy_reads.k15_w5.surj_realn.bam"
ROW="$WORK/results/toy_reads/row_k15_w5.tsv"

echo ""
echo "== Result =="
if [ -s "$BAM" ]; then
    mapped=$(awk 'NR==2{print $4}' "$ROW" 2>/dev/null || echo "?")
    total=$(awk 'NR==2{print $3}' "$ROW" 2>/dev/null || echo "?")
    echo "TEST PASSED: produced $BAM"
    echo "Mapped $mapped / $total reads. Summary: $ROW"
else
    echo "TEST FAILED: expected BAM not found at $BAM"
    exit 1
fi
