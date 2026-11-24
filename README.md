# pangenome-mapping

Map short-read sequencing data to the **HPRC v1.1 pangenome** with `vg giraffe`
and produce a BAM on linear GRCh38 coordinates (giraffe surjects on the fly,
`-o BAM`) plus mapping statistics. No indel realignment step is applied.

The pipeline uses **personalized haplotype sampling** (`vg haplotypes`): for each
sample it builds a smaller, sample-specific graph from k-mer content before
mapping, which improves both speed and accuracy over mapping to the full graph.

## Pipeline overview

```
FASTQ (one or more runs per sample, single- or paired-end)
   │
   ├─ 1. KMC              count k-mers in all runs of the sample (k=29)
   ├─ 2. vg haplotypes    sample ONE personalized graph (.gbz) per sample
   │                      (k=29, w=11; 4 candidate haplotypes by default)
   ├─ 3. vg index         distance index
   ├─ 4. vg minimizer     minimizer index (mapping k/w)
   ├─ 5. vg giraffe -o BAM  map each run, output BAM on GRCh38 paths
   │     rename_bam_stream.py   strip the "GRCh38#0#" contig prefix
   │     samtools sort; per-run BAMs merged per sample
   └─ 6. samtools stats   mapping-quality summary (row_k<k>_w<w>.tsv, MAPQ dist)
```

Output for a sample lands in `results/<sample>/`, ending in
`<sample>.k<k>_w<w>.bam` (+ `.bai`) and the stats TSVs.

## Requirements

These tools must be on your `PATH` at run time (versions the pipeline was
developed with in parentheses):

| tool | used for |
|------|----------|
| [`vg`](https://github.com/vgteam/vg) (≥ 1.71) | haplotype sampling, indexing, giraffe |
| [`kmc`](https://github.com/refresh-bio/KMC) | k-mer counting |
| [`samtools`](https://www.htslib.org/) | sort / index / stats |
| `python3` | `rename_bam_stream.py` |

## Setup

```bash
git clone https://github.com/admixVIE/pangenome-mapping.git
cd pangenome-mapping

# 1. Download the large public reference files (~25 GB, one time).
#    Needs samtools on PATH to build the .fai/.dict.
bash references/download_references.sh

# 2. Create your local config from the template and edit the paths/environment.
cp config.example.sh config.sh
$EDITOR config.sh
```

## Configuration

All machine-specific settings live in `config.sh` (git-ignored, so your paths
are never committed). Copy the template and edit it:

```bash
cp config.example.sh config.sh
$EDITOR config.sh
```

Set two things (plus optional tuning, below):

- **Paths** — `PANGENOME_ROOT` (your `references/` dir) and `REF_FA`.
- **Environment** — the `load_environment()` function, so that `vg`, `kmc`,
  `samtools` and `python3` are on the `PATH`. On a
  cluster with modules use e.g. `ml load KMC vg SAMtools`;
  with conda use `conda activate <your-env>`. `config.example.sh` has examples
  for both.

> Outside the original LISC cluster, creating `config.sh` is required — the
> built-in fallbacks use LISC-specific module names.

## Usage

Input is always a **runs table** (TSV with header `sample  run  fq1  fq2`; leave
`fq2` empty for single-end). A sample may have several runs: each run is mapped
separately and the BAMs are merged per sample.

```bash
bash mapping/map_to_pangenome.sh -S Bon002 --runs runs.tsv -t 8

# override the mapping minimizer or the number of sampled haplotypes
bash mapping/map_to_pangenome.sh -S Bon002 --runs runs.tsv -t 8 -k 15 -w 5 --num-haplotypes 32
```

### Mapping k/w and haplotypes

- **Mapping minimizer k/w:** set `MAP_KMER` / `MAP_WINDOW` in `config.sh`, or pass
  `-k` / `-w`. If neither is set, **vg giraffe's defaults (k=29, w=11)** are used.
  The k/w used is part of the output name (`<sample>.k29_w11.bam`).
- **Haplotype sampling** always uses k=29, w=11 (fixed) and **4 candidate
  haplotypes** by default (`NUM_HAPLOTYPES` in `config.sh` or `--num-haplotypes`).
  Changing it changes the personalized graph, so a sample must be re-run from
  scratch (delete its `results/<sample>/`).

Run with no arguments to print the full option list (`-r`, `-ref`, `-g`, `-h`,
`-p`, `-pre`, ...). On a SLURM cluster, edit and submit
[`mapping/run_mapping.slurm`](mapping/run_mapping.slurm):

```bash
sbatch mapping/run_mapping.slurm
```

The script is **resumable**: each step is skipped if its output already exists,
so a re-run continues where a failed run stopped.

## Quick test

A tiny synthetic pangenome ships with the repo so you can run the whole pipeline
in seconds without downloading the large references — handy for checking your
environment is set up correctly:

```bash
bash mapping/test/test_run.sh
# -> TEST PASSED: ... toy_reads.k15_w5.bam  (Mapped 4000 / 4000 reads)
```

See [mapping/test/README.md](mapping/test/README.md) for details.

## Simulated reads and benchmarking

`reads_simulation/` holds the gargammel-based pipeline that generates damaged,
short single-end read sets (3 read-length distributions x 3 damage levels x
4 coverages) from the HG002 assembly, and `benchmarking/` the hap.py scripts
used to score linear vs pangenome mapping on them. Both need their paths set
via environment variables / a config file; see the README in each folder.

## Repository layout

```
pangenome-mapping/
├── README.md
├── config.example.sh          # copy to config.sh (git-ignored) and edit
├── mapping/
│   ├── map_to_pangenome.sh     # the pipeline
│   ├── rename_bam_stream.py    # strips the graph's contig-name prefix
│   ├── run_mapping.slurm       # SLURM submission wrapper
│   └── test/                   # end-to-end smoke test + tiny toy pangenome
├── reads_simulation/           # gargammel + ART simulation of damaged short reads (see its README)
├── benchmarking/               # hap.py scoring + summary/plot scripts (see its README)
└── references/
    ├── download_references.sh  # fetches hg38 + HPRC graph (not in git)
    └── hprc-v1.1-mc-grch38.paths_list.txt   # small, ships with the repo
```

## Outputs

For each run in `results/<sample>/`:

- `<sample>.k<k>_w<w>.bam` (+ `.bai`) — final alignments on GRCh38 coordinates
- `*_samtools.stats`, `*_flagstat.txt` — full samtools statistics
- `row_k<k>_w<w>.tsv` — one-line summary: total / mapped / unmapped / MAPQ≥30
- `*_mapq_distribution.tsv` — MAPQ histogram for plotting

## References

- Pangenome graph: HPRC v1.1 Minigraph-Cactus, GRCh38 reference
  ([human-pangenomics](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes%2Ffreeze%2Ffreeze1%2Fminigraph-cactus%2Fhprc-v1.1-mc-grch38%2F))
- Haplotype sampling method: Sirén et al., *Personalized pangenome references*
  ([bioRxiv 2023.12.13.571553](https://www.biorxiv.org/content/10.1101/2023.12.13.571553))
- `vg`: [vgteam/vg](https://github.com/vgteam/vg)
