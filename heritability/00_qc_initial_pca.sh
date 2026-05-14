#!/bin/bash
# 00_qc_initial_pca.sh — Initial QC and PCA on the pre-QC'd WGS cohort
#
# Input:
#   - Pre-QC'd PLINK fileset from the upstream qc/ pipeline
#     (LCR filter, genotype QC, call rate, KING relatedness, sex check, HWE)
#
# Output:
#   - <prefix>_qc.{bed,bim,fam}        post-QC PLINK fileset
#   - PPS_pca100.eigenvec / .eigenval  100 PCs for case-control matching
#
# Notes:
#   - MAF > 1% threshold matches paper Methods (common-variant heritability).
#   - LD pruning (--indep-pairwise 50 5 0.1) is required before PCA so
#     correlated SNPs don't dominate the eigendecomposition.

set -euo pipefail

BUCKET="${WORKSPACE_BUCKET}"
PREFIX="PPS_431cases_10418controls_lcr_filter_genotypeQC_callrate_king_sex_hwe_wgs_plink"

# Copy input PLINK files
for ext in bed bim fam; do
  gsutil -m cp "${BUCKET}/data/${PREFIX}.${ext}" .
done

# Variant + genotype QC
plink --noweb --allow-no-sex --memory 80000 \
  --bfile "${PREFIX}" \
  --maf 0.01 --geno 0.98 \
  --make-bed --out "${PREFIX}_qc"

# LD pruning
plink --noweb --allow-no-sex --memory 80000 \
  --bfile "${PREFIX}_qc" \
  --indep-pairwise 50 5 0.1 \
  --out "${PREFIX}_qc"

# 100 PCs on pruned set (used for matching in step 02; only PC1-5 are actually used)
plink --noweb --allow-no-sex --memory 90000 --threads 200 \
  --bfile "${PREFIX}_qc" \
  --extract "${PREFIX}_qc.prune.in" \
  --pca 100 \
  --out PPS_pca100

echo "100 PCs written to PPS_pca100.eigenvec"
