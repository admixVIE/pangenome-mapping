# pangenome-mapping

Map short-read sequencing data to the **HPRC v1.1 pangenome** with `vg giraffe`,
project the alignments back onto the linear GRCh38 reference (surjection), and
produce an indel-realigned, analysis-ready BAM plus mapping statistics.

The pipeline uses **personalized haplotype sampling** (`vg haplotypes`): for each
sample it builds a smaller, sample-specific graph from k-mer content before
mapping, which improves both speed and accuracy over mapping to the full graph.

## Pipeline overview

```
FASTQ (single- or paired-end)
   │
   ├─ 1. KMC              count k-mers in the reads
   ├─ 2. vg haplotypes    sample a personalized graph (.gbz) from the HPRC graph
   ├─ 3. vg index         distance index
   ├─ 4. vg minimizer     minimizer index
   ├─ 5. vg giraffe       map reads to the personalized graph  → GAF
   ├─ 6. vg surject       project GAF onto GRCh38 paths        → BAM
   │     rename_bam_stream.py   strip the "GRCh38#0#" contig prefix
   │     bamleftalign + samtools sort
   ├─ 8. GATK3 RealignerTargetCreator   find indel-realignment targets
   ├─ 9. ABRA2                          indel realignment       → final BAM
   └─  samtools stats     mapping-quality summary (row_k<k>_w<w>.tsv, MAPQ dist)
```

Output for a sample lands in `results/<sample>_k<k>_w<w>/`, ending in
`*.surj_realn.bam` (+ `.bai`) and the stats TSVs.

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
git clone <your-gitlab-url>/pangenome-mapping.git
cd pangenome-mapping

# 1. Download the large public reference files (~30 GB, one time).
#    Needs samtools on PATH to build the .fai/.dict.
bash references/download_references.sh

# 2. Create your local config from the template and edit the paths/environment.
cp config.example.sh config.sh
$EDITOR config.sh
```

`config.sh` holds everything machine-specific (reference location, the GATK 3.8
jar, and a `load_environment()` function that loads your modules / conda envs).
It is git-ignored. If you run on the original LISC cluster the built-in fallbacks
already work, so `config.sh` is optional there.

## Usage

```bash
# single-end
bash scripts/map_to_pangenome.sh -s1 reads.fq.gz -t 8

# paired-end
bash scripts/map_to_pangenome.sh -s1 reads_R1.fq.gz -s2 reads_R2.fq.gz -t 8

# change vg k-mer / window parameters
bash scripts/map_to_pangenome.sh -s1 reads.fq.gz -t 8 -k 15 -w 5
```

Run with no arguments to print the full option list (`-r`, `-ref`, `-g`, `-h`,
`-p`, `-pre`, ...). On a SLURM cluster, edit and submit
[`scripts/run_mapping.slurm`](scripts/run_mapping.slurm):

```bash
sbatch scripts/run_mapping.slurm
```

The script is **resumable**: each step is skipped if its output already exists,
so a re-run continues where a failed run stopped.

## Repository layout

```
pangenome-mapping/
├── README.md
├── config.example.sh          # copy to config.sh (git-ignored) and edit
├── scripts/
│   ├── map_to_pangenome.sh     # the pipeline
│   ├── rename_bam_stream.py    # strips the graph's contig-name prefix
│   └── run_mapping.slurm       # SLURM submission wrapper
├── references/
│   ├── download_references.sh  # fetches hg38 + HPRC graph + ABRA2 (not in git)
│   └── hprc-v1.1-mc-grch38.paths_list.txt   # small, ships with the repo
└── test/                       # small runnable example (see test/README.md)
```

## Outputs

For each run in `results/<sample>_k<k>_w<w>/`:

- `*.surj_realn.bam` (+ `.bai`) — final indel-realigned alignments on GRCh38
- `*_samtools.stats`, `*_flagstat.txt` — full samtools statistics
- `row_k<k>_w<w>.tsv` — one-line summary: total / mapped / unmapped / MAPQ≥30
- `*_mapq_distribution.tsv` — MAPQ histogram for plotting

## References

- Pangenome graph: HPRC v1.1 Minigraph-Cactus, GRCh38 reference
  ([human-pangenomics](https://s3-us-west-2.amazonaws.com/human-pangenomics/index.html?prefix=pangenomes%2Ffreeze%2Ffreeze1%2Fminigraph-cactus%2Fhprc-v1.1-mc-grch38%2F))
- Haplotype sampling method: Sirén et al., *Personalized pangenome references*
  ([bioRxiv 2023.12.13.571553](https://www.biorxiv.org/content/10.1101/2023.12.13.571553))
- `vg`: [vgteam/vg](https://github.com/vgteam/vg)
