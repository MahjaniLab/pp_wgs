"""
02_filter_flagged.py — Remove AoU-flagged samples

The All of Us QC team publishes a flagged_samples.tsv listing samples that
failed central QC metrics (variant-count outliers, contamination, etc.).

Input:
  - ${WORKSPACE_BUCKET}/data/cohort.mt (from step 01)
  - AoU flagged_samples.tsv (env var AOU_FLAGGED_SAMPLES_PATH)
    For v7: gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/qc/flagged_samples.tsv

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_flagged_passed.mt

"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
FLAGGED_PATH = os.getenv("AOU_FLAGGED_SAMPLES_PATH")

IN_MT  = f"{BUCKET}/data/cohort.mt"
OUT_MT = f"{BUCKET}/data/cohort_flagged_passed.mt"

mt = hl.read_matrix_table(IN_MT)

flagged = hl.import_table(FLAGGED_PATH).key_by("s")
mt = mt.annotate_cols(**flagged[mt.s])

# Drop samples with any QC flag
before = mt.count_cols()
mt = mt.filter_cols(hl.len(mt.qc_metrics_filters) == 0)
after = mt.count_cols()

print(f"Removed {before - after} flagged samples")
print(f"Samples remaining: {after}")

mt.write(OUT_MT, overwrite=True)
