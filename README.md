# pangenome-mapping

Map short-read sequencing data to the **HPRC v1.1 pangenome** with `vg giraffe`,
project the alignments back onto the linear GRCh38 reference (surjection), and
produce an indel-realigned, analysis-ready BAM plus mapping statistics.

The pipeline uses **personalized haplotype sampling** (`vg haplotypes`): for each
sample it builds a smaller, sample-specific graph from k-mer content before
mapping, which improves both speed and accuracy over mapping to the full graph.

## Pipeline overview

```
FASTQ (one or more runs per sample, single- or paired-end)
   │
   ├─ 1. KMC              count k-mers in all runs of the sample (k=29)
   ├─ 2. vg haplotypes    sample ONE personalized graph (.gbz) per sample (k=29, w=11)
   ├─ 3. vg index         distance index
   ├─ 4. vg minimizer     minimizer index (mapping k/w, default 15/5)
   ├─ 5. vg giraffe       map each run to the personalized graph → one GAF per run
   ├─ 6. vg surject       project each GAF onto GRCh38 paths     → one BAM per run, merged per sample
   │     rename_bam_stream.py   strip the "GRCh38#0#" contig prefix
   │     bamleftalign + samtools sort
   ├─ 8. GATK3 RealignerTargetCreator   find indel-realignment targets
   ├─ 9. ABRA2                          indel realignment       → final BAM
   └─  samtools stats     mapping-quality summary (row_k<k>_w<w>.tsv, MAPQ dist)
```

Output for a sample lands in `results/<sample>/`, ending in
`<sample>.k<k>_w<w>.surj_realn.bam` (+ `.bai`) and the stats TSVs.

## Requirements

These tools must be on your `PATH` at run time (versions the pipeline was
developed with in parentheses):

| tool | used for |
|------|----------|
| [`vg`](https://github.com/vgteam/vg) (≥ 1.71) | haplotype sampling, indexing, giraffe, surject |
| [`kmc`](https://github.com/refresh-bio/KMC) | k-mer counting |
| [`samtools`](https://www.htslib.org/) | sort / index / stats |
| [`bedtools`](https://bedtools.readthedocs.io/) | target BED construction |
| `bamleftalign` (from [freebayes](https://github.com/freebayes/freebayes)) | left-align indels |
| **GATK 3.8** + **Java 8** | `RealignerTargetCreator` |
| [ABRA2](https://github.com/mozack/abra2) 2.23 (jar) | indel realignment |

> GATK 3.8 is required specifically — `RealignerTargetCreator` was removed in
> GATK 4. It needs a Java 8 runtime.

## Setup

```bash
git clone https://github.com/admixVIE/pangenome-mapping.git
cd pangenome-mapping

# 1. Download the large public reference files (~30 GB, one time).
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

Set two things:

- **Paths** — `PANGENOME_ROOT` (your `references/` dir), and `GATK3_JAR` / `JAVA8`.
- **Environment** — the `load_environment()` function, so that `vg`, `kmc`,
  `samtools`, `bedtools`, `bamleftalign` and `java` are on the `PATH`. On a
  cluster with modules use e.g. `ml load KMC vg SAMtools BEDTools freebayes`;
  with conda use `conda activate <your-env>`. `config.example.sh` has examples
  for both.

> Outside the original LISC cluster, creating `config.sh` is required — the
> built-in fallbacks use LISC-specific module names.

## Usage

```bash
# single-end
bash mapping/map_to_pangenome.sh -s1 reads.fq.gz -t 8

# paired-end
bash mapping/map_to_pangenome.sh -s1 reads_R1.fq.gz -s2 reads_R2.fq.gz -t 8

# several runs per sample (runs.tsv: sample, run, fq1, fq2; fq2 empty for single-end)
bash mapping/map_to_pangenome.sh -S Bon002 --runs runs.tsv -t 8

# change the mapping minimizer k / w (default 15 / 5; haplotype sampling is fixed at 29 / 11)
bash mapping/map_to_pangenome.sh -S Bon002 --runs runs.tsv -t 8 -k 29 -w 11
```

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
# -> TEST PASSED: ... toy_reads.k15_w5.surj_realn.bam  (Mapped 3993 / 4000 reads)
```

See [mapping/test/README.md](mapping/test/README.md) for details.

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
└── references/
    ├── download_references.sh  # fetches hg38 + HPRC graph + ABRA2 (not in git)
    └── hprc-v1.1-mc-grch38.paths_list.txt   # small, ships with the repo
```

## Outputs

For each run in `results/<sample>/`:

- `<sample>.k<k>_w<w>.surj_realn.bam` (+ `.bai`) — final indel-realigned alignments on GRCh38
- `*_samtools.stats`, `*_flagstat.txt` — full samtools statistics
- `row_k<k>_w<w>.tsv` — one-line summary: total / mapped / unmapped / MAPQ≥30
- `*_mapq_distribution.tsv` — MAPQ histogram for plotting

## References

- Pangenome graph: HPRC v1.1 Minigraph-Cactus, GRCh38 reference
  ([human-pangenomics](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes%2Ffreeze%2Ffreeze1%2Fminigraph-cactus%2Fhprc-v1.1-mc-grch38%2F))
- Haplotype sampling method: Sirén et al., *Personalized pangenome references*
  ([bioRxiv 2023.12.13.571553](https://www.biorxiv.org/content/10.1101/2023.12.13.571553))
- `vg`: [vgteam/vg](https://github.com/vgteam/vg)
