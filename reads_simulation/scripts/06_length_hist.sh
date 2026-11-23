#!/bin/bash
set -euo pipefail
cd "$(dirname "$0")/.."
for c in short mid long; do
  ( zcat reads/final/$c/high/cov1/merged.fq.gz | awk 'NR%4==2{n[length($0)]++} END{for(l in n) print l"\t"n[l]}' | sort -n > figures/read_length_counts_$c.tsv ) &
done
wait
