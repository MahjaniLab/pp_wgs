"""
03_filter_lcr.py — Remove variants in low-complexity regions

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_flagged_passed.mt (from step 02)
  - LCR-hs38.bed at ${WORKSPACE_BUCKET}/data/LCR-hs38.bed
    Heng Li's LCR resource for GRCh38; pre-load to the bucket.

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_lcr.mt

Reference:
  Li H. Toward better understanding of artifacts in variant calling from
  high-coverage samples. Bioinformatics 2014.
"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
LCR_BED = f"{BUCKET}/data/LCR-hs38.bed"
IN_MT   = f"{BUCKET}/data/cohort_flagged_passed.mt"
OUT_MT  = f"{BUCKET}/data/cohort_lcr.mt"

mt = hl.read_matrix_table(IN_MT)

lcr = hl.import_bed(LCR_BED, reference_genome="GRCh38")
mt = mt.filter_rows(hl.is_defined(lcr[mt.locus]), keep=False)

print(f"Variants after LCR removal: {mt.count_rows()}")

mt.write(OUT_MT, overwrite=True)
