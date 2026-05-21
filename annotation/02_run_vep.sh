#!/bin/bash
# 02_run_vep.sh — Functional annotation with Ensembl VEP + LOFTEE
#
# Runs VEP in offline cache mode with the LOFTEE and Conservation plugins to
# annotate consequence, canonical transcript, gene symbol, SIFT/PolyPhen, and
# high-confidence loss-of-function (LoF) calls.
#
# Run on a host with the VEP cache and LOFTEE resources installed (local
# workstation or HPC node), not on the AoU workbench. Download vep_input.txt
# from the bucket first. Input is VCF-style (chrom pos . ref alt).
#
# Input:
#   - vep_input.txt (from step 01)
# Output:
#   - vep_output.txt
#
# --- Set these to your local VEP / LOFTEE installation ---

VEP=/path/to/ensembl-vep/vep

CACHE_DIR=/path/to/ensembl-vep/cache

PLUGIN_DIR=/path/to/loftee

FASTA=/path/to/hg38.fa

HUMAN_ANCESTOR=/path/to/human_ancestor.fa.gz

CONSERVATION_SQL=/path/to/loftee.sql

GERP_BW=/path/to/gerp_conservation_scores.homo_sapiens.GRCh38.bw

IN=vep_input.txt

OUT=vep_output.txt

"${VEP}" \
  -i "${IN}" \
  -o "${OUT}" \
  --offline --cache --dir_cache "${CACHE_DIR}" \
  --assembly GRCh38 \
  --fasta "${FASTA}" \
  --force_overwrite \
  --no_stats \
  --everything \
  --variant_class \
  --allele_number \
  --pick_allele \
  --minimal \
  --canonical \
  --ccds \
  --biotype \
  --symbol \
  --sift b \
  --polyphen b \
  --dir_plugin "${PLUGIN_DIR}" \
  --plugin LoF,loftee_path:"${PLUGIN_DIR}",human_ancestor_fa:"${HUMAN_ANCESTOR}",conservation_file:"${CONSERVATION_SQL}",gerp_bigwig:"${GERP_BW}",gerp_database:"${GERP_BW}" \
  --plugin Conservation,"${GERP_BW}"

echo "VEP annotation written to ${OUT}"
