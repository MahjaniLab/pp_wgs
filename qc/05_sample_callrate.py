"""
05_sample_callrate.py — Sample call rate filter

Removes samples with call rate below (mean - 3 SD). This adaptive threshold
accommodates platform- and dataset-specific call rate distributions.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_genotypeQC.mt (from step 04)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_callrate.mt
  - ${WORKSPACE_BUCKET}/data/sample_callrate.tsv

Expected exclusions: ~43 samples (paper Methods, v7).
"""

import os
import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT  = f"{BUCKET}/data/cohort_genotypeQC.mt"
OUT_MT = f"{BUCKET}/data/cohort_callrate.mt"
CR_TSV = f"{BUCKET}/data/sample_callrate.tsv"

mt = hl.read_matrix_table(IN_MT)
mt = hl.sample_qc(mt)

mt.cols().select("sample_qc").export(CR_TSV)

stats = mt.aggregate_cols(hl.agg.stats(mt.sample_qc.call_rate))
threshold = stats.mean - 3 * stats.stdev

print(f"Mean call rate: {stats.mean:.4f}")
print(f"SD:             {stats.stdev:.4f}")
print(f"Threshold (mean - 3SD): {threshold:.4f}")

before = mt.count_cols()
mt = mt.filter_cols(mt.sample_qc.call_rate >= threshold)
after = mt.count_cols()

print(f"Removed {before - after} samples; {after} remain")

mt.write(OUT_MT, overwrite=True)
