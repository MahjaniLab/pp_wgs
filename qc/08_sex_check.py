"""
08_sex_check.py — chrX sex check

The cohort is restricted to parous women, so all samples should genotype as
female. Samples with F-statistic > 0.4 likely represent sex-chromosome
aneuploidies or mislabeled samples and are removed.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_unrelated.mt (from step 07)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_sex_pass.mt
  - ${WORKSPACE_BUCKET}/data/sex_check_results.tsv

"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT  = f"{BUCKET}/data/cohort_unrelated.mt"
OUT_MT = f"{BUCKET}/data/cohort_sex_pass.mt"
RESULTS_TSV = f"{BUCKET}/data/sex_check_results.tsv"

F_THRESHOLD = 0.4   # F > 0.4 -> fail (cohort is all female)

mt = hl.read_matrix_table(IN_MT)

# Restrict to chrX, common high-quality variants
chrX = [hl.parse_locus_interval("chrX", reference_genome="GRCh38")]
mt_x = hl.filter_intervals(mt, chrX)
mt_x = hl.variant_qc(mt_x)
mt_x = mt_x.filter_rows(
    (mt_x.variant_qc.AF[1] >= 0.1) & (mt_x.variant_qc.call_rate >= 0.98)
)

imputed = hl.impute_sex(mt_x.GT)
imputed.export(RESULTS_TSV)

# Identify samples failing
fail_ids = imputed.filter(imputed.f_stat > F_THRESHOLD).s.collect()
print(f"Sex check failures (F > {F_THRESHOLD}): {len(fail_ids)}")

fail_set = hl.literal(set(fail_ids))
mt = mt.filter_cols(~fail_set.contains(mt.s))

print(f"Samples after sex check: {mt.count_cols()}")

mt.write(OUT_MT, overwrite=True)
