#!/bin/bash
# 05_run_pcgc.sh — PCGC heritability estimation with LDAK6
#
# Two-step procedure for each GRM (autosomes, chrX):
#   1. --adjust-grm: regress out covariates from the GRM
#   2. --pcgc:       PCGC regression with disease prevalence
#
# Input:
#   - HumDef.grm.*, HumDef_xchr.grm.*   (from step 04)
#   - pca20.eigenvec                    (from step 03; 20 PCs)
#   - pheno.txt                         (from step 03)
#
# Output:
#   - pcgc_autosome.pcgc
#   - pcgc_xchr.pcgc
#   - pcgc_summary.txt
#
# Parameters:
#   prevalence = 0.0015
#     Lifetime prevalence of postpartum psychosis among parous women,
#     from the Swedish national register cohort (n=1,648,759); Table S2.
#
# Software:
#   LDAK6: https://dougspeed.com/ldak/
#   PCGC reference: Golan D, Lander ES, Rosset S. PNAS 2014.
#   https://doi.org/10.1073/pnas.1419031111
#
# Paper-reported values (Methods, Table 1):
#   h2_autosomal = 0.434 (SE 0.224)
#   h2_X         = 0.022 (SE 0.055)
#   h2_total     = 0.456 (SE 0.222)

set -euo pipefail

PREV=0.0015
COVAR=pca20.eigenvec   # 20 PCs per paper Methods

# ----- Autosomal -----
./ldak6.linux --adjust-grm HumDef.covar --grm HumDef --covar "${COVAR}"

./ldak6.linux \
  --pcgc pcgc_autosome \
  --pheno pheno.txt \
  --prevalence "${PREV}" \
  --grm HumDef.covar \
  --covar "${COVAR}"

# ----- chrX -----
./ldak6.linux --adjust-grm HumDef_xchr.covar --grm HumDef_xchr --covar "${COVAR}"

./ldak6.linux \
  --pcgc pcgc_xchr \
  --pheno pheno.txt \
  --prevalence "${PREV}" \
  --grm HumDef_xchr.covar \
  --covar "${COVAR}"

# ----- Summary -----
{
  echo "=== PCGC heritability summary ==="
  echo "Prevalence: ${PREV}"
  echo "Covariates: ${COVAR} (20 PCs)"
  echo ""
  echo "--- Autosomal ---"
  cat pcgc_autosome.pcgc
  echo ""
  echo "--- chrX ---"
  cat pcgc_xchr.pcgc
  echo ""
  echo "Paper-reported (Methods, Table 1):"
  echo "  h2_autosomal = 0.434 (SE 0.224)"
  echo "  h2_X         = 0.022 (SE 0.055)"
  echo "  h2_total     = 0.456 (SE 0.222)"
} > pcgc_summary.txt

cat pcgc_summary.txt
