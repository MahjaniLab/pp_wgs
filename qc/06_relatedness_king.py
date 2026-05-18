"""
06_relatedness_king.py — Sample relatedness via Hail KING

Computes kinship coefficients with Hail's KING implementation, exports related
pairs (kinship > 0.0884, equivalent to 2nd-degree or closer), and automatically
drops samples that appear in many related pairs (>15 connections).

The remaining pairs require manual review — see README for the curation
protocol. Once the user has prepared manual_drop.tsv, run 07_apply_manual_drop.py.

Input:
  - ${WORKSPACE_BUCKET}/data/cohort_callrate.mt (from step 05)

Output:
  - ${WORKSPACE_BUCKET}/data/related_pairs.tsv          (all pairs above threshold)
  - ${WORKSPACE_BUCKET}/data/high_connection_drop.tsv   (auto-dropped samples)
  - ${WORKSPACE_BUCKET}/data/cohort_unrelated_auto.mt   (MT after auto-drop)

Threshold rationale:
  KING > 0.0884 -> 2nd-degree relatives or closer.
  Equivalent to the conventional PLINK PI_HAT > 0.2 cutoff.

Reference:
  Manichaikul A et al. Bioinformatics 2010 (KING-robust).
"""

import os
import hail as hl
import pandas as pd

hl.init(default_reference="GRCh38")

BUCKET = os.getenv("WORKSPACE_BUCKET")
IN_MT  = f"{BUCKET}/data/cohort_callrate.mt"
OUT_MT = f"{BUCKET}/data/cohort_unrelated_auto.mt"
PAIRS_TSV = f"{BUCKET}/data/related_pairs.tsv"
DROP_TSV  = f"{BUCKET}/data/high_connection_drop.tsv"

KING_THRESHOLD = 0.0884   # 2nd-degree or closer
MAX_CONNECTIONS = 15      # auto-drop samples with this many related partners

mt = hl.read_matrix_table(IN_MT)

# Restrict to common, high-quality autosomal biallelic SNPs for KING
mt_qc = hl.variant_qc(mt)
mt_king = mt_qc.filter_rows(
    (mt_qc.variant_qc.call_rate >= 0.98) &
    (mt_qc.variant_qc.AF[1] >= 0.01) &
    (hl.len(mt_qc.alleles) == 2) &
    (mt_qc.locus.in_autosome())
)

# LD prune
pruned = hl.ld_prune(mt_king.GT, r2=0.1, bp_window_size=500_000)
mt_king = mt_king.filter_rows(hl.is_defined(pruned[mt_king.row_key]))

# KING kinship MT (samples x samples)
king = hl.king(mt_king.GT)

# Extract upper-triangle pairs above threshold
pairs = king.entries()
pairs = pairs.filter((pairs.s_1 < pairs.s) & (pairs.phi > KING_THRESHOLD))
pairs = pairs.rename({"s_1": "IID1", "s": "IID2", "phi": "kinship"})
pairs.select("IID1", "IID2", "kinship").export(PAIRS_TSV)

print(f"Related pairs (KING > {KING_THRESHOLD}): {pairs.count()}")

# Count connections per sample
pairs_df = pd.read_csv(PAIRS_TSV, sep="\t")
counts = pd.concat([pairs_df["IID1"], pairs_df["IID2"]]).value_counts()
high_conn = counts[counts > MAX_CONNECTIONS].index.tolist()

print(f"Samples with >{MAX_CONNECTIONS} related partners: {len(high_conn)}")
pd.DataFrame({"sample_id": high_conn}).to_csv(DROP_TSV, sep="\t", index=False)

# Auto-drop high-connection samples
drop_set = hl.literal(set(high_conn))
mt = mt.filter_cols(~drop_set.contains(mt.s))

print(f"Samples after auto-drop: {mt.count_cols()}")
mt.write(OUT_MT, overwrite=True)

print(f"""
==============================================================
NEXT STEP: manual curation
  1. Open {PAIRS_TSV}.
  2. Exclude rows where either IID is already in {DROP_TSV}.
  3. For each remaining pair, choose ONE sample to drop:
       - prefer cases over controls
       - within same phenotype, prefer the higher call rate sample
  4. Save chosen drops to {BUCKET}/data/manual_drop.tsv
     with a sample_id column header.
  5. Run 07_apply_manual_drop.py.
==============================================================
""")
