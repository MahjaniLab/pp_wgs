# 02_match_optmatch.R — 1:10 case-control matching with optmatch
#
# Matches each case to 10 controls on PC1-5 using Euclidean distance.
# Matched cohort is intersected with the NFE filter in step 03.
#
# Input:
#   - <prefix>_qc.fam      (phenotype in column 6: PLINK 2=case, 1=control)
#   - PPS_pca100.eigenvec  (100 PCs from step 00)
#
# Output:
#   - PPS_WGS_matched_1to10_result.txt: matched cohort with PCs + match strata
#   - matched_PC1_PC2.png: QC plot
#
# Reference:
#   Hansen BB, Klopfer SO. JCGS 2006. https://doi.org/10.1198/106186006X137047

suppressPackageStartupMessages({
  library(optmatch)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)
options(optmatch_max_problem_size = Inf)

PREFIX <- "PPS_431cases_10418controls_lcr_filter_genotypeQC_callrate_king_sex_hwe_wgs_plink_qc"

fam <- read.table(paste0(PREFIX, ".fam"), header = FALSE)
colnames(fam) <- c("FID", "IID", "Father", "Mother", "sex", "type")

pcs <- read.table("PPS_pca100.eigenvec", header = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:100))

input <- merge(fam, pcs[, c("IID", paste0("PC", 1:5))], by = "IID", sort = FALSE)

# PLINK pheno coding (2=case, 1=control) -> (1=case, 0=control) for matching
input$pheno <- as.integer(ifelse(input$type == 2, 1, 0))

cat("Before matching:\n")
cat("  Cases:   ", sum(input$pheno == 1), "\n")
cat("  Controls:", sum(input$pheno == 0), "\n\n")

# Euclidean distance on PC1-5
distances <- match_on(
  pheno ~ PC1 + PC2 + PC3 + PC4 + PC5,
  data = input, method = "euclidean"
)

# 1:10 matching
match10 <- pairmatch(distances, data = input, controls = 10)

matched <- cbind(input, matches = match10)
matched <- matched[!is.na(matched$matches), ]

cat("After 1:10 matching:\n")
print(table(matched$pheno))

write.table(matched, "PPS_WGS_matched_1to10_result.txt",
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")

# QC plot: verify case/control overlap in PC space
matched$pheno_lab <- ifelse(matched$pheno == 1, "case", "control")

ggplot(matched, aes(PC1, PC2, color = pheno_lab, shape = pheno_lab)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c(case = "red", control = "black")) +
  theme_classic() +
  ggtitle("Matched cohort: PC1 vs PC2")

ggsave("matched_PC1_PC2.png", width = 6, height = 5, dpi = 150)
