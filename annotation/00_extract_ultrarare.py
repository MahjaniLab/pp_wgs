"""
00_extract_ultrarare.py — Extract ultra-rare variants

Ultra-rare variants are defined as allele count <= 5 within the cohort AND
<= 5 (or absent) in the gnomAD v3.1 non-neuro subset. These variants are
carried forward into functional annotation and rare-variant analyses.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_matched.mt (from matching/)
  - gnomAD v3.1.2 sites HT (public)

Output:
  - ${WORKSPACE_BUCKET}/data/cohort_ultrarare.mt
"""

import os

import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")

IN_MT  = f"{BUCKET}/data/cohort_matched.mt"

OUT_MT = f"{BUCKET}/data/cohort_ultrarare.mt"

GNOMAD_HT = "gs://gcp-public-data--gnomad/release/3.1.2/ht/genomes/gnomad.genomes.v3.1.2.sites.ht"

mt = hl.read_matrix_table(IN_MT)

# Cohort allele count 1-5 (ultra-rare within cohort)

mt = hl.variant_qc(mt)

mt = mt.filter_rows((mt.variant_qc.AC[1] > 0) & (mt.variant_qc.AC[1] <= 5))

# gnomAD non-neuro AC <= 5 (or absent in gnomAD)

gnomad = hl.read_table(GNOMAD_HT)

non_neuro_idx = hl.eval(gnomad.freq_index_dict["non_neuro-adj"])

mt = mt.annotate_rows(gnomad_non_neuro_AC=gnomad[mt.row_key].freq[non_neuro_idx].AC)

mt = mt.filter_rows(

    hl.is_missing(mt.gnomad_non_neuro_AC) | (mt.gnomad_non_neuro_AC <= 5)

)

print(f"Ultra-rare variants: {mt.count_rows():,}")

mt.write(OUT_MT, overwrite=True)
