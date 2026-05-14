#!/bin/bash
# 01_filter_nfe.sh — Restrict to Non-Finnish European (NFE) ancestry
#
# Uses ancestry predictions from a PCA projection onto 1000 Genomes
# super-populations. Samples with NFE probability > 0.8 are retained.
# Finnish-ancestry samples are excluded due to bottleneck-driven relatedness
# that would inflate the GRM.
#
# Input:
#   - PPS_368cases_9238controls_1000G_merge_pca_predicted_populations.txt
#     Columns: sample_id  ...  predicted_population  probability
#
# Output:
#   - NFE.txt: tab-separated FID IID for plink --keep

set -euo pipefail

BUCKET="${WORKSPACE_BUCKET}"
ANCESTRY_FILE="PPS_368cases_9238controls_1000G_merge_pca_predicted_populations.txt"

gsutil -m cp "${BUCKET}/data/${ANCESTRY_FILE}" .

awk '$4 == "NFE" && $5 > 0.8 {print $1"\t"$1}' "${ANCESTRY_FILE}" > NFE.txt

echo "NFE samples: $(wc -l < NFE.txt)"
