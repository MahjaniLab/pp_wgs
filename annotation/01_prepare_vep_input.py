"""
01_prepare_vep_input.py — Export ultra-rare variants to VCF-style VEP input

Output format (tab-separated, no header, no 'chr' prefix):
  chrom  pos  .  ref  alt

e.g.  1  69360  .  C  T

VEP auto-detects this as VCF-style input.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_ultrarare.mt (from step 00)

Output:
  - ${WORKSPACE_BUCKET}/data/vep_input.txt
    (download to the VEP host before running step 02)
"""

import os

import hail as hl

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")

IN_MT = f"{BUCKET}/data/cohort_ultrarare.mt"

OUT   = f"{BUCKET}/data/vep_input.txt"

mt = hl.read_matrix_table(IN_MT)

rows = mt.rows()

rows = rows.annotate(

    chrom=rows.locus.contig.replace("chr", ""),   # chr1 -> 1

    pos=rows.locus.position,

    vid=hl.str("."),

    ref=rows.alleles[0],

    alt=rows.alleles[1],

)

rows = rows.key_by().select("chrom", "pos", "vid", "ref", "alt")

rows.export(OUT, header=False)

print(f"VEP input written to {OUT}")
