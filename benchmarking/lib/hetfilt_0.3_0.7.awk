# Keep a heterozygous call only if alt/(ref+alt) read support is in [0.3, 0.7].
# Input: single-sample, biallelic VCF stream. Homozygous calls and records without GT/AD pass unchanged.
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
