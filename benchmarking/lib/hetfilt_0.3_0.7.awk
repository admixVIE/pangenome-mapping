# hetfilt_0.3_0.7.awk -- reconstructed 2026-09-16 (original deleted during a
# scratch-space cleanup, no backup existed). Called on a single-sample,
# already-biallelic (post `bcftools norm -m-`) VCF stream as the last stage of
# the "confidence" set: raw AND FILTER==PASS AND INFO/DP>=5 AND het 0.3-0.7.
#
# Header lines and any record without a usable GT/AD pass through untouched.
# Homozygous genotypes (0/0, 1/1, ...) pass through unconditionally -- this
# filter only targets heterozygous allele-balance skew. A het call is kept
# only if AD[alt allele]/(AD[ref allele]+AD[alt allele]) falls in [0.3, 0.7];
# outside that range it's dropped as a likely miscall (e.g. near-homozygous
# support masquerading as het).
BEGIN { FS = OFS = "\t" }
/^#/ { print; next }
{
  ngt = split($9, fmt, ":")
  gtIdx = 0; adIdx = 0
  for (i = 1; i <= ngt; i++) {
    if (fmt[i] == "GT") gtIdx = i
    if (fmt[i] == "AD") adIdx = i
  }
  if (gtIdx == 0) { print; next }
  split($10, smp, ":")
  gt = smp[gtIdx]
  gsub(/\|/, "/", gt)
  n = split(gt, alleles, "/")
  if (n < 2 || alleles[1] == "." || alleles[2] == ".") { print; next }
  if (alleles[1] == alleles[2]) { print; next }
  if (adIdx == 0) { print; next }
  split(smp[adIdx], ad, ",")
  a1 = alleles[1] + 0; a2 = alleles[2] + 0
  d1 = ad[a1 + 1] + 0; d2 = ad[a2 + 1] + 0
  tot = d1 + d2
  if (tot == 0) { print; next }
  ratio = d2 / tot
  if (ratio < 0.3 || ratio > 0.7) next
  print
}
