# Test data

A tiny, fully synthetic pangenome (~260 KB, committed to the repo) for running
the pipeline end-to-end in a few seconds — no ~30 GB reference download needed.

## Run it

```bash
bash test/test_run.sh
```

This runs the full `map_to_pangenome.sh` pipeline on the toy data and checks that
a non-empty indel-realigned BAM with mapped reads is produced. It uses the same
tools as a real run (vg, kmc, samtools, bedtools, bamleftalign, GATK 3.8 + Java 8,
ABRA2), so it also confirms your `config.sh` / environment is set up correctly.
Set that up first (see the main [README](../README.md)); on the original LISC
cluster the built-in fallbacks are used.

Expected output (numbers may vary slightly by tool version):

```
TEST PASSED: produced .../toy_reads.k15_w5.surj_realn.bam
Mapped 3993 / 4000 reads.
```

Outputs are written to `test/work/` (git-ignored).

## What's in `data/`

| file | description |
|------|-------------|
| `chr1.fa` (+ `.fai`, `.dict`) | 40 kb synthetic reference contig `chr1` |
| `toy.giraffe.gbz` | small pangenome graph (3 samples, 19 SNPs) |
| `toy.hapl` | haplotype-sampling index |
| `toy.paths_list.txt` | reference path name (`chr1`) |
| `toy_reads.fq.gz` | 4000 single-end 100 bp reads simulated from the graph |

All of it is randomly generated — there is **no real human sequence** here.

The toy graph and reference both use the plain contig name `chr1`, so no name
prefix needs stripping (the test passes `-pre "none#"`, a no-op). On real HPRC
data the prefix is `GRCh38#0#`.

## Regenerating the data

The fixtures are committed, so you don't need to. To rebuild them from scratch
(requires vg, samtools, HTSlib, python3):

```bash
bash test/make_test_data.sh
```
