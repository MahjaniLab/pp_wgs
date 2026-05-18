"""
01_extract_cohort_mt.py — Subset AoU WGS MatrixTable to the PPS cohort

Annotates each sample with phenotype (1=case, 0=control).

Input:
  - AoU srWGS MatrixTable (env var WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH)
  - ${WORKSPACE_BUCKET}/data/cohort_phenotype.csv (from step 00)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort.mt
"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
WGS_MT_PATH = os.getenv("WGS_ACAF_THRESHOLD_MULTI_HAIL_PATH")
PHENO_FILE = f"{BUCKET}/data/cohort_phenotype.csv"
OUT_MT = f"{BUCKET}/data/cohort.mt"

mt = hl.read_matrix_table(WGS_MT_PATH)

pheno = hl.import_table(
    PHENO_FILE, delimiter=",", types={"phenotype": hl.tint32}
).key_by("person_id")

mt = mt.filter_cols(hl.is_defined(pheno[mt.s]))
mt = mt.annotate_cols(phenotype=pheno[mt.s].phenotype)

print(f"Cases:    {mt.filter_cols(mt.phenotype == 1).count_cols()}")
print(f"Controls: {mt.filter_cols(mt.phenotype == 0).count_cols()}")

mt.write(OUT_MT, overwrite=True)
print(f"Wrote {OUT_MT}")
