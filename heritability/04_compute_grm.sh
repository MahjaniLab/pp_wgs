#!/bin/bash
# 04_compute_grm.sh — Compute LDAK-Thin kinship matrices
#
# Uses LDAK6 with --power -.25 (LDAK-Thin model), which downweights
# common SNPs and gives less biased h2 estimates than the GCTA model
# for traits with rare-variant contributions.
#
# Input:
#   - PPS_wgs_matched_qc_autosome_prune.{bed,bim,fam}
#   - PPS_wgs_matched_qc_xchr_prune.{bed,bim,fam}
#
# Output:
#   - HumDef.grm.{bin,id,N}        (autosomal)
#   - HumDef_xchr.grm.{bin,id,N}   (chrX)
#
# Software:
#   LDAK6 by Doug Speed: https://dougspeed.com/ldak/
#   Reference: Speed D et al. Nat Genet 2017; 2020.

set -euo pipefail

BUCKET="${WORKSPACE_BUCKET}"

# Fetch LDAK6 binary if not already present
if [ ! -x ./ldak6.linux ]; then
  gsutil -m cp "${BUCKET}/data/ldak6.linux" .
  chmod +x ldak6.linux
fi

# Autosomal GRM (LDAK-Thin)
./ldak6.linux \
  --calc-kins-direct HumDef \
  --bfile PPS_wgs_matched_qc_autosome_prune \
  --power -.25 \
  --max-threads 100

# chrX GRM
./ldak6.linux \
  --calc-kins-direct HumDef_xchr \
  --bfile PPS_wgs_matched_qc_xchr_prune \
  --power -.25 \
  --max-threads 100

echo "GRMs written: HumDef (autosomal), HumDef_xchr (chrX)"
