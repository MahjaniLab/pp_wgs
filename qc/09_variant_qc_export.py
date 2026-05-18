"""
09_variant_qc_export.py — Final variant QC and PLINK export

Final variant-level QC and conversion to PLINK format for downstream pipelines
(heritability, matching, burden testing).

  - Variant call rate >= 0.9
  - HWE p-value >= 1e-12
  - PLINK export with phenotype encoded as 2 (case) / 1 (control)

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_sex_pass.mt (from step 08)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_final.mt
  - ${WORKSPACE_BUCKET}/data/cohort_final.{bed,bim,fam}

"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT     = f"{BUCKET}/data/cohort_sex_pass.mt"
OUT_MT    = f"{BUCKET}/data/cohort_final.mt"
OUT_PLINK = f"{BUCKET}/data/cohort_final"

mt = hl.read_matrix_table(IN_MT)

mt = hl.variant_qc(mt)
mt = mt.filter_rows(mt.variant_qc.call_rate >= 0.9)
mt = mt.filter_rows(mt.variant_qc.p_value_hwe >= 1e-12)

n_var  = mt.count_rows()
n_samp = mt.count_cols()
n_case = mt.filter_cols(mt.phenotype == 1).count_cols()
n_ctrl = mt.filter_cols(mt.phenotype == 0).count_cols()

print(f"Final variants: {n_var:,}")
print(f"Final samples:  {n_samp} ({n_case} cases, {n_ctrl} controls)")

mt.write(OUT_MT, overwrite=True)

# PLINK export (phenotype: 2=case, 1=control)
mt = mt.annotate_cols(plink_pheno=hl.if_else(mt.phenotype == 1, 2, 1))
hl.export_plink(mt, OUT_PLINK, ind_id=mt.s, pheno=mt.plink_pheno)

print(f"PLINK files written to {OUT_PLINK}.{{bed,bim,fam}}")
