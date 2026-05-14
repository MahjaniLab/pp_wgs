#!/bin/bash
# 03_finalize_matched.sh — Build final analysis cohort for LDAK
#
# Intersects matched samples (step 02) with NFE ancestry (step 01), then:
#   - subsets the PLINK fileset
#   - applies MAF > 1% filter (paper Methods)
#   - splits autosomes (1-22) and chrX (PLINK chr 23)
#   - LD-prunes each separately
#   - re-computes 20 PCs on the matched, pruned autosomal data (PCGC covariates)
#   - extracts pheno.txt for LDAK
#
# Input:
#   - <prefix>_qc.{bed,bim,fam}              (from step 00)
#   - PPS_WGS_matched_1to10_result.txt       (from step 02)
#   - NFE.txt                                (from step 01)
#
# Output:
#   - PPS_wgs_matched_qc_autosome_prune.{bed,bim,fam}
#   - PPS_wgs_matched_qc_xchr_prune.{bed,bim,fam}
#   - pca20.eigenvec                         (20 PCs for PCGC covariates)
#   - pheno.txt                              (FID IID phenotype)

set -euo pipefail

PREFIX_IN="PPS_431cases_10418controls_lcr_filter_genotypeQC_callrate_king_sex_hwe_wgs_plink_qc"

# 1. matched ∩ NFE
awk 'NR>1 {print $2"\t"$2}' PPS_WGS_matched_1to10_result.txt | sort -u > matched_ids.txt
sort -u NFE.txt > NFE_sorted.txt
comm -12 matched_ids.txt NFE_sorted.txt > matched_nfe.txt
echo "Matched ∩ NFE: $(wc -l < matched_nfe.txt)"

# 2. Subset + final QC (MAF > 1% per paper Methods)
plink --noweb --allow-no-sex --memory 90000 --threads 200 \
  --bfile "${PREFIX_IN}" \
  --keep matched_nfe.txt \
  --maf 0.01 --geno 0.98 \
  --make-bed --out PPS_wgs_matched_qc

# 3. Split autosomes / chrX (PLINK chr 23 = X)
plink --noweb --allow-no-sex --memory 90000 \
  --bfile PPS_wgs_matched_qc --autosome \
  --make-bed --out PPS_wgs_matched_qc_autosome

plink --noweb --allow-no-sex --memory 90000 \
  --bfile PPS_wgs_matched_qc --chr 23 \
  --make-bed --out PPS_wgs_matched_qc_xchr

# 4. LD-prune each
plink --noweb --allow-no-sex --memory 90000 \
  --bfile PPS_wgs_matched_qc_autosome \
  --indep-pairwise 50 5 0.1 --out prune_autosome

plink --noweb --allow-no-sex \
  --bfile PPS_wgs_matched_qc_xchr \
  --indep-pairwise 50 5 0.1 --out prune_xchr

plink --noweb --allow-no-sex \
  --bfile PPS_wgs_matched_qc_autosome \
  --extract prune_autosome.prune.in \
  --make-bed --out PPS_wgs_matched_qc_autosome_prune

plink --noweb --allow-no-sex \
  --bfile PPS_wgs_matched_qc_xchr \
  --extract prune_xchr.prune.in \
  --make-bed --out PPS_wgs_matched_qc_xchr_prune

# 5. 20 PCs on matched, pruned autosomal data (PCGC covariates per paper Methods)
plink --noweb --allow-no-sex --memory 90000 \
  --bfile PPS_wgs_matched_qc_autosome \
  --extract prune_autosome.prune.in \
  --pca 20 --out pca20

# 6. Pheno file for LDAK (FID IID phenotype; LDAK reads PLINK 2/1 coding)
awk '{print $1"\t"$2"\t"$6}' PPS_wgs_matched_qc_autosome_prune.fam > pheno.txt

echo "Final cohort: $(wc -l < PPS_wgs_matched_qc_autosome_prune.fam) samples"
echo "Autosomal variants: $(wc -l < PPS_wgs_matched_qc_autosome_prune.bim)"
echo "chrX variants:      $(wc -l < PPS_wgs_matched_qc_xchr_prune.bim)"
