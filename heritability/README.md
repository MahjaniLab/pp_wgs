# WGS-based heritability estimation

We use PCGC regression (Golan et al., PNAS 2014) as implemented in LDAK6 (Speed et al., Nat Genet 2017; 2020) to estimate SNP heritability of postpartum psychosis. PCGC is designed for binary case–control traits with strong ascertainment — controls were oversampled in our matched cohort, and postpartum psychosis has a population prevalence of only ~0.15%.

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 0 | 00_qc_initial_pca.sh | PLINK QC + 100 PCs for matching |
| 1 | 01_filter_nfe.sh | NFE ancestry filter (1000G PCA projection, NFE prob > 0.8) |
| 2 | 02_match_optmatch.R | 1:10 case-control matching on PC1–5 (R optmatch) |
| 3 | 03_finalize_matched.sh | Subset matched ∩ NFE, autosome/chrX split, LD prune, 20 PCs |
| 4 | 04_compute_grm.sh | LDAK-Thin GRM --calc-kins-direct --power -.25) |
| 5 | 05_run_pcgc.sh | LDAK --adjust-grm + --pcgc |

## Key parameters

- Cohort: 198 EUR (NFE) cases + 2,013 EUR controls (combined AoU v7 + v8)
- Variants: MAF > 1%, autosomes + chrX, LD-pruned (common-variant heritability)
- GRM model: LDAK-Thin --power -.25) — downweights common SNPs
- Covariates: first 20 PCs (computed after QC and ancestry restriction)
- Prevalence: 0.0015 (Swedish register estimate among parous women; Table S2)

## Why these choices?

- LDAK-Thin over GCTA: less biased h2 estimates for traits with rare-variant contributions (Speed et al. 2017).
- NFE-only: PCGC is sensitive to population stratification on the GRM. Finnish-ancestry samples are excluded due to bottleneck-driven relatedness.
- 20 PCs: captures fine-scale within-EUR structure without overcorrecting.
- Autosomes vs. chrX separately: PCGC treats them differently; results are combined post-hoc (autosomal 43.4% + chrX 2.2% = 45.6% total).

## Software

- LDAK6: https://dougspeed.com/ldak/
- PLINK 1.9: https://www.cog-genomics.org/plink/
- R packages: optmatch, dplyr, ggplot2

## References

- Golan D, Lander ES, Rosset S. PNAS 2014. https://doi.org/10.1073/pnas.1419064111
- Speed D et al. Nat Genet 2017. https://doi.org/10.1038/ng.3865
- Speed D et al. Nat Genet 2020. https://doi.org/10.1038/s41588-020-0600-y
- Hansen BB, Klopfer SO. JCGS 2006 optmatch).

## Reproducing the paper's result

Run scripts 00–05 in order. Expected output from 05_run_pcgc.sh:

```
h2_autosomal: 0.434 (SE 0.224)
h2_X:         0.022 (SE 0.055)
h2_total:     0.456 (SE 0.222)
```
