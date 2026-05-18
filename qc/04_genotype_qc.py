"""
04_genotype_qc.py — Genotype-level QC

Entry-level filters for high-confidence genotype calls:
  - FT == "PASS" (or missing)
  - GQ >= 25
  - Allele balance:
      hom_ref   : AB <= 0.1
      het       : AB >= 0.3 AND pab >= 1e-9 (two-sided binomial test)
      hom_var   : AB >= 0.9
  - Variants with zero called genotypes after filtering are removed.

Notes:
  AB = AD[alt] / sum(AD); pab is the binomial test of AD[alt] vs 0.5.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_lcr.mt (from step 03)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_genotypeQC.mt
"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT  = f"{BUCKET}/data/cohort_lcr.mt"
OUT_MT = f"{BUCKET}/data/cohort_genotypeQC.mt"

mt = hl.read_matrix_table(IN_MT)

# FT field
mt = mt.filter_entries((mt.FT == "PASS") | (hl.is_missing(mt.FT)))

# AB + GQ + pab
ab  = mt.AD[1] / hl.sum(mt.AD)
pab = hl.binom_test(mt.AD[1], hl.sum(mt.AD), 0.5, "two-sided")

keep = (
    ((mt.GT.is_hom_ref()) & (ab <= 0.1) & (mt.GQ >= 25)) |
    ((mt.GT.is_hom_ref()) & hl.is_missing(mt.AD) & (mt.GQ >= 25)) |
    ((mt.GT.is_het())     & (ab >= 0.3) & (mt.GQ >= 25) & (pab >= 1e-9)) |
    ((mt.GT.is_hom_var()) & (ab >= 0.9) & (mt.GQ >= 25))
)

mt = mt.filter_entries(keep)

# Drop sites with no called genotypes
mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.n_called > 0)

print(f"Variants after genotype QC: {mt.count_rows()}")

mt.write(OUT_MT, overwrite=True)
