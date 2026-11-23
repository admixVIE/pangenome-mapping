#!/bin/bash
# ==============================================================================
# make_test_data.sh — regenerate the toy pangenome fixtures in test/data/.
# ------------------------------------------------------------------------------
# You normally do NOT need to run this: the fixtures are committed to the repo.
# It is kept for transparency/reproducibility. Requires vg, samtools, HTSlib
# (bgzip/tabix) and python3 on PATH.
#
# What it builds (all tiny, fully synthetic — no real human data):
#   * chr1.fa (+ .fai, .dict)   40 kb random reference contig "chr1"
#   * a phased VCF with 19 SNPs across 3 samples -> a small pangenome graph
#   * toy.giraffe.gbz           graph (vg autoindex, giraffe workflow)
#   * toy.hapl                  haplotype-sampling index (vg haplotypes -H)
#   * toy.paths_list.txt        reference path name ("chr1")
#   * toy_reads.fq.gz           4000 single-end 100 bp reads (vg sim)
# ==============================================================================
set -euo pipefail
OUT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/data"
mkdir -p "$OUT"; cd "$OUT"

python3 - <<'PY'
import random
random.seed(42)
seq=''.join(random.choice('ACGT') for _ in range(40000))
with open('chr1.fa','w') as f:
    f.write('>chr1\n')
    for i in range(0,len(seq),60): f.write(seq[i:i+60]+'\n')
alt={'A':'G','G':'A','C':'T','T':'C'}
samples=['s1','s2','s3']; positions=list(range(2000,40000,2000))[:19]
with open('vars.vcf','w') as f:
    f.write('##fileformat=VCFv4.2\n##FILTER=<ID=PASS,Description="pass">\n')
    f.write('##contig=<ID=chr1,length=40000>\n')
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+'\t'.join(samples)+'\n')
    for j,pos in enumerate(positions):
        ref=seq[pos-1]; a=alt[ref]
        gts=[f'{(j+si)%2}|{1 if (j+si+1)%3==0 else 0}' for si in range(len(samples))]
        f.write(f'chr1\t{pos}\t.\t{ref}\t{a}\t60\tPASS\t.\tGT\t'+'\t'.join(gts)+'\n')
PY

bgzip -f vars.vcf && tabix -p vcf vars.vcf.gz
samtools faidx chr1.fa
samtools dict chr1.fa -o chr1.dict

vg autoindex --workflow giraffe -r chr1.fa -v vars.vcf.gz -p toy -t 4
vg gbwt -r toy.ri -Z toy.giraffe.gbz
vg haplotypes -t 4 --subchain-length 1000 -d toy.dist -r toy.ri -H toy.hapl toy.giraffe.gbz
vg sim -x toy.giraffe.gbz -n 4000 -l 100 -e 0.005 -a | vg view -X -a - | gzip > toy_reads.fq.gz
printf 'chr1\n' > toy.paths_list.txt

# keep only the files the test needs
rm -f vars.vcf.gz vars.vcf.gz.tbi graphref.fa.fai \
      toy.dist toy.ri toy.shortread.withzip.min toy.shortread.zipcodes
echo "Done. Fixtures written to $OUT"
