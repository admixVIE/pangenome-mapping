# Benchmarking (hap.py)

Scores variant calls made from the simulated reads (see `../reads_simulation/`)
against the HG002 T2T-Q100 v1.1 small-variant truth set with **hap.py**
(vcfeval engine, run through an apptainer image).

| script | what it does |
|---|---|
| `scripts/happy_array.sbatch` | hap.py on `confidence.vcf.gz` of one finished set (SLURM array; line N of `$LIST` = `<arm> <set>`) |
| `scripts/happy_dp1_array.sbatch` | 1x only: rebuild confidence from `raw.vcf.gz` with DP>=1, then hap.py |
| `scripts/make_summary_table.py` | collect `summary_happy.tsv` into the pangenome-vs-linear table |
| `scripts/plot_summary_figures.py` | precision / recall / F1 figures |
| `scripts/plot_read_length_kde.py`, `plot_mapdamage_combined.py` | read-length and damage-profile figures |
| `lib/` | `lib_harmonize.sh` (contig-name / ALT=* harmonisation), het-AB awk filter, `rtg-wrapped-mem` (caps rtg's Java heap) |

## Setup

```bash
cp config.example.sh config.sh && $EDITOR config.sh   # truth dir, hap.py image, rtg-tools, BENCH_ROOT
cd benchmarking
sbatch --array=0-N --export=ALL,LIST=/path/to/list.txt scripts/happy_array.sbatch
```

The python scripts read `BENCH_ROOT` from the environment (the results working
directory). Not tested outside the original cluster.
