# Shared harmonization functions -- sourced by 04_check_format_compat.sh and
# the {with_wave,without_wave,k29w11}_02_bench_truvari.sh /
# _03_bench_happy.sh scripts so the fix lives in exactly one place.
#
# Findings this encodes (see README.md "Format-compatibility findings"):
#   - our VCF: PanSN CHROM (GRCh38#0#chr1) -> strip to bare chr1
#   - ALT=* placeholders (157,310 of them in the truth) must be dropped with
#     -e 'ALT="*"' (GIAB's own README's -e 'ALT="."' is wrong for this file --
#     verified empirically)
#   - restrict both sides to the 24 primary contigs the truth header declares
#
# CORRECTED 2026-09-09: this used to also pre-filter the truth to FILTER=="."
# here, reasoning that truvari's --passonly does a literal string match on
# "PASS" and would zero out a truth file whose confident marker is "." not
# "PASS". That reasoning was wrong -- it was based on truvari's --help text,
# not its actual code. Checked variant_record.py's is_filtered(): with no
# explicit values, it treats an EMPTY filter set (which is what "." parses to
# in pysam/htslib -- confirmed empirically against this exact truth file) as
# NOT filtered, same as literal PASS. So --passonly already does the right
# thing on its own, on both truth and comp; the manual pre-filter here was
# redundant. Removed so 05_truvari_bench.sh's --passonly matches GIAB's
# documented command exactly, rather than duplicating the same logic twice.

BCFTOOLS=${BCFTOOLS:-bcftools}
TRUVARI=${TRUVARI:-$HOME/.local/bin/truvari}
PANSN_PREFIX="GRCh38#0#"

# harmonize_truth <raw_truth.vcf.gz> <out_prefix>
# Writes <out_prefix>.vcf.gz(.tbi): ALT!="*" only, sorted+indexed. Confidence
# filtering (FILTER=="."  or PASS) is left to truvari's own --passonly.
harmonize_truth() {
  local raw="$1" outp="$2"
  [ -s "${outp}.vcf.gz" ] && return 0
  "$BCFTOOLS" view -e 'ALT="*"' "$raw" -Oz -o "${outp}.vcf.gz"
  "$BCFTOOLS" index -t "${outp}.vcf.gz"
}

# harmonize_query <raw_query.vcf.gz> <out_prefix> <chrom_map.tsv>
# chrom_map.tsv: 2 columns, "GRCh38#0#chr1<TAB>chr1" etc (bcftools annotate --rename-chrs format)
# Strips PanSN prefix, drops ALT=*, restricts to primary contigs, sorts+indexes.
harmonize_query() {
  local raw="$1" outp="$2" map="$3" primary_contigs="$4"
  # -T (targets-file, streams through) not -R (regions-file, needs an index on
  # the actual file) -- -R on piped/unindexed BCF fails with "Could not
  # retrieve index file for '-'" (confirmed 2026-09-06 against the real pilot
  # VCF). Also found empirically on that same run: this pipeline's real
  # output is ALREADY bare `chr1`, not PanSN `GRCh38#0#chr1` -- the vcfwave/
  # norm chain in 09_normalize_vcf.sh evidently strips it before we ever see
  # it, contrary to that script's own "STILL OPEN" comment. The rename step
  # below is therefore a harmless no-op on current output, kept in case that
  # changes.
  "$BCFTOOLS" annotate --rename-chrs "$map" "$raw" -Ou \
    | "$BCFTOOLS" view -e 'ALT="*"' -Ou \
    | "$BCFTOOLS" view -T "$primary_contigs" -Ou \
    | "$BCFTOOLS" sort -Oz -o "${outp}.vcf.gz"
  "$BCFTOOLS" index -t "${outp}.vcf.gz"
  fix_inv_flag_header "${outp}.vcf.gz"
}

# fix_inv_flag_header <vcf.gz>
# In-place fix, only if needed: this pipeline's INV annotation writes
# INFO/INV=YES but declares it Type=Flag,Number=0 in the header -- a Flag
# field must never carry "=value" per the VCF spec. Harmless while nothing
# reads the value back out, but truvari refine's final GA4GH-consolidation
# step round-trips every record through pysam, which enforces the spec and
# crashes with "ValueError: Flag values must be: True, False, None, 1, 0" on
# the first INV=YES record it hits (confirmed 2026-09-11 on cond B raw/het
# and cond C raw -- these carry an inversion that survives into the FN/FP
# set; cond A has no INV field at all, cond B/C dp5 happened not to have an
# INV record land in a refined region). Same class of bug as the AD
# Number=. issue -- fixed the same way, via reheader (not `annotate -h`,
# which does not let you change a field's declared Type).
fix_inv_flag_header() {
  local vcf="$1"
  "$BCFTOOLS" view -h "$vcf" | grep -q '^##INFO=<ID=INV,Number=0,Type=Flag' || return 0
  local hdr; hdr=$(mktemp)
  "$BCFTOOLS" view -h "$vcf" \
    | sed 's/^##INFO=<ID=INV,Number=0,Type=Flag,\(.*\)>$/##INFO=<ID=INV,Number=1,Type=String,\1>/' \
    > "$hdr"
  "$BCFTOOLS" reheader -h "$hdr" -o "${vcf}.reheader" "$vcf"
  mv "${vcf}.reheader" "$vcf"
  rm -f "$hdr"
  "$BCFTOOLS" index -tf "$vcf"
}

# build_chrom_map <surject_paths.txt> <out.tsv>
# surject_paths.txt has one PanSN path per line, e.g. "GRCh38#0#chr1"
build_chrom_map() {
  local paths="$1" out="$2"
  # Guard added 2026-09-18: `> "$out"` truncates the destination BEFORE awk
  # even opens "$paths", so a missing/moved source silently wipes out an
  # existing, working chrom_map.tsv instead of just failing loudly.
  [ -s "$paths" ] || { echo "ERROR: build_chrom_map: source not found: $paths (leaving $out untouched)" >&2; return 1; }
  awk -v pfx="$PANSN_PREFIX" '{orig=$0; sub("^"pfx,"",$0); print orig"\t"$0}' "$paths" > "$out"
}

# primary_contigs_bed <truth.vcf.gz> <out.bed>
# One whole-contig row per ##contig declared in the truth header (the 24
# primary chroms) -- used to restrict the query VCF to the same contig set.
primary_contigs_bed() {
  local vcf="$1" out="$2" fai="$3"
  "$BCFTOOLS" view -h "$vcf" | grep '^##contig' | sed -E 's/.*ID=([^,>]+).*/\1/' > "${out}.names"
  awk 'NR==FNR{want[$1]=1; next} want[$1]{print $1"\t0\t"$2}' "${out}.names" "$fai" > "$out"
}
