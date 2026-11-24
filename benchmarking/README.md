# Benchmarking (hap.py)

Scores variant calls made from the simulated reads (see `../reads_simulation/`)
against the HG002 T2T-Q100 v1.1 small-variant truth set with **hap.py**
(vcfeval engine, run through an apptainer image).

| script | what it does |
|---|---|
| `scripts/happy_array.sbatch` | one SLURM array script: builds the confidence call set from `raw.vcf.gz` (PASS, `DP >= MIN_DP`, het allele balance 0.3-0.7) and scores it with hap.py; `MIN_DP` defaults to 5 |
| `lib/` | `hetfilt_0.3_0.7.awk` (het allele-balance filter), `rtg-wrapped-mem` (caps rtg's Java heap; optional, drop `--engine-vcfeval-path` if you have enough RAM) |

Summary tables and plots are not included: `postproc/results/summary_happy.tsv`
is a plain TSV (arm, set, coverage, threshold tag, hap.py columns) to plot as you like.

## Setup

```bash
cd benchmarking
cp config.example.sh config.sh && $EDITOR config.sh   # truth dir, hap.py image, rtg-tools, BENCH_ROOT
sbatch --array=0-N --export=ALL,LIST=/path/to/list.txt scripts/happy_array.sbatch
# other depth threshold, e.g. 1x sets:  --export=ALL,LIST=...,MIN_DP=1
```

Not tested outside the original cluster.
