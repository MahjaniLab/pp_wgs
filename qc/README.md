# Sample and variant QC

QC pipeline for the All of Us short-read WGS cohort, producing the post-QC PLINK fileset that feeds into downstream analyses (heritability, matching, burden testing).

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 0 | 00_define_cohort.py | Build PPS case/control list from AoU SQL extracts |
| 1 | 01_extract_cohort_mt.py | Subset AoU srWGS MT to cohort, annotate phenotype |
| 2 | 02_filter_flagged.py | Drop AoU-flagged samples (variant count outliers, contamination, etc.) |
| 3 | 03_filter_lcr.py | Remove variants in low-complexity regions |
| 4 | 04_genotype_qc.py | Genotype-level FT / AB / GQ / pab filters |
| 5 | 05_sample_callrate.py | Drop samples with call rate < mean − 3SD |
| 6 | 06_relatedness_king.py | Hail KING; auto-drop samples with >15 related partners |
| — | (manual review) | Curate manual_drop.tsv from remaining related pairs |
| 7 | 07_apply_manual_drop.py | Apply manual drop list |
| 8 | 08_sex_check.py | chrX F-stat sex check (cohort is all female) |
| 9 | 09_variant_qc_export.py | Variant call rate + HWE + PLINK export |

## Cohort definition

Cases: women with mania or psychosis (ICD-9 / ICD-10 codes per Table S1 of the paper) within 0–3 months of their first live birth. Controls: parous women without any psychiatric diagnosis. The combined v7 cohort yielded 431 cases and 10,418 candidate controls before QC.

## Sample QC summary (paper Methods, v7)

| Filter | Samples removed |
|---|---|
| AoU-flagged | 31 |
| Call rate | 43 |
| Relatedness (KING) | 189 |
| Sex mismatch | 7 |
| Remaining | 305 cases + 8,736 controls |

## Manual relatedness review

06_relatedness_king.py exports two files: related_pairs.tsv (all pairs with KING > 0.0884) and high_connection_drop.tsv (samples with > 15 related partners, auto-dropped). For pairs still remaining after the auto-drop, manual review is required:

1. Open related_pairs.tsv and remove rows where either sample appears in high_connection_drop.tsv.
2. For each remaining pair, pick one sample to drop. Preference order:
   - Keep cases over controls.
   - Within the same phenotype, keep the higher call rate sample.
3. Save the final drop list as manual_drop.tsv with a sample_id column header.
4. Run 07_apply_manual_drop.py.

## Key thresholds

- KING kinship: 0.0884 (2nd-degree or closer; equivalent to PLINK PI_HAT > 0.2)
- Genotype QC: GQ ≥ 25; allele balance hom_ref ≤ 0.1 / het 0.3–0.7 (binom pab ≥ 1e-9) / hom_var ≥ 0.9
- Sample call rate: mean − 3 standard deviations (adaptive)
- Variant call rate: ≥ 0.9
- HWE: p ≥ 1e-12
- Sex F-stat: F > 0.4 → fail (cohort is all female)

## Software

- Hail (≥ 0.2.130): https://hail.is/
- PLINK 1.9 (downstream only)

## References

- Manichaikul A et al. Bioinformatics 2010 (KING).
- Karczewski KJ et al. Nature 2020 (gnomAD QC).
- Li H. Bioinformatics 2014 (LCR-hs38.bed).
