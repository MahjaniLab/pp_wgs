"""
07_apply_manual_drop.py — Apply manual relatedness drop list

After reviewing related_pairs.tsv (step 06) and curating manual_drop.tsv,
remove the chosen samples to produce the final unrelated cohort.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_unrelated_auto.mt (from step 06)
  - ${WORKSPACE_BUCKET}/data/manual_drop.tsv          (user-curated; column: sample_id)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_unrelated.mt
"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT     = f"{BUCKET}/data/cohort_unrelated_auto.mt"
DROP_FILE = f"{BUCKET}/data/manual_drop.tsv"
OUT_MT    = f"{BUCKET}/data/cohort_unrelated.mt"

mt   = hl.read_matrix_table(IN_MT)
drop = hl.import_table(DROP_FILE).key_by("sample_id")

before = mt.count_cols()
mt = mt.filter_cols(~hl.is_defined(drop[mt.s]))
after = mt.count_cols()

print(f"Removed {before - after} samples via manual review")
print(f"Samples remaining: {after}")

mt.write(OUT_MT, overwrite=True)
