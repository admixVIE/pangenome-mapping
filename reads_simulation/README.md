# Read simulation (experiment 1: linear vs pangenome mapping, HG002)

**Reference script (commands and parameters of all conditions in one place):** `scripts/simulate_all_conditions_REFERENCE.sh`
The pipeline was actually run as SLURM steps `01a_fragments`, `01_simulate`, `02_trim_merge`, `03_subsample` (+ `04_art_split/part/join` for five conditions); parameters in `scripts/config.sh`.

## Design
36 single-end read sets = 3 read-length distributions x 3 damage levels x 4 coverages (1x, 5x, 7x, 10x), all simulated from the HG002 T2T-Q100 v1.1 assembly (maternal + paternal haplotypes + chrM, 6.0 Gb). One simulation per condition (no replicates, fixed seeds). Only reads that were adapter-trimmed and merged are used for mapping.

## Steps
1. **Fragments** (gargammel `fragSim`). Fragments are drawn from the assembly in proportion to contig length (assembly given as one haploid "endo" file), length 30-291 bp.
   - short: empirical distribution from a real library, median 46 / mean 55 bp
   - mid: empirical distribution from a real library, median 67 / mean 84 bp
   - long: lognormal (meanlog 4.6, sdlog 0.4) truncated to 30-291 bp, median ~100 / mean 107 bp
   - Number of fragments = 10 x 2,999,712,359 bp / mean length: short 545.6 M, mid 356.2 M, long 280.2 M. One fragment set per length distribution, shared by the three damage levels (so damage levels are paired).
2. **Damage** (`deamSim`, Briggs model, double-stranded library, `-damage v,l,d,s`). Parameters calibrated so that the C->T rate at the first 5' base is about 7 / 15 / 25 %:
   - low `0.024,0.36,0.0031,0.206`, mid `0.024,0.36,0.0066,0.447`, high `0.024,0.36,0.0108,0.76`
   - Check on the final data with mapDamage (10x, mean of the three length classes): C->T at 5' position 1 = 0.068 / 0.144 / 0.241, G->A at 3' position 1 = 0.067 / 0.142 / 0.237 (`summary/mapdamage/mapdamage_combined.png`).
3. **Adapters** (`adptSim`, Illumina TruSeq-type adapters), so short fragments read through into adapter.
4. **Sequencing** (ART Illumina `art_illumina`, profile HS25, paired-end 2 x 125 bp, `-amp -na -p -c 1 -qs 0 -qs2 0`). Default ART error model: substitutions and small insertions/deletions are included, so some INDEL false positives can come from simulated sequencing errors. For five conditions (short/low, short/mid, mid/high, long/low, long/mid) ART was run in 20 parts with different seeds because single runs stalled; parameters identical.
5. **Trimming and merging** (`leeHom --ancientdna`, version 1.2.18). Read pairs are trimmed and merged; only fragments up to about 239 bp (2 x 125 bp minus the minimum overlap) can be merged, so the longest fragments (240-291 bp) are lost. Merged reads: 94-96 % of fragments. Unmerged pairs exist but are not used.
6. **Coverage** (`seqtk sample -s 321`). Reads at 10x are randomly sub-sampled with fractions 0.7, 0.5, 0.1 (same seed: the 1x set is nested in 5x, in 7x, in 10x).

## Resulting read sets (merged reads)
| length | fragments simulated (10x) | merged reads at 5x | median / mean merged length (bp) | 5 / 95 % quantile |
|---|---|---|---|---|
| short | 545.6 M | 262.9 M | 46 / 54.5 | 31 / 109 |
| mid | 356.2 M | 168.0 M | 66 / 80.6 | 34 / 180 |
| long | 280.2 M | 132.9 M | 99 / 105.2 | 51 / 184 |
(merged lengths from the first 20 M reads of each low-damage 5x set; `summary/read_lengths/`)

## Caveats
- "Coverage" is nominal (fragments x mean length / 3.0 Gb). After mapping and filtering (MQ, duplicates) the depth at called sites is lower (median HaplotypeCaller DP about 0.8x of nominal at 5x; about 0.7x at 10x).
- Autosomes receive twice the depth of X/Y per contig pair, as in a male.
- One replicate per condition; no confidence intervals from repeated simulations.
- Read names carry the true origin (`contig:strand:start:end`), but the assembly coordinates are not those of hg38; mapping accuracy per read was not used here.
