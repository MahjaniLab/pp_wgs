# Genetic Architecture of Postpartum Psychosis 

Code and analysis pipelines accompanying:

> Jung S, Caballero M, Kępińska A, Smout S, Munk-Olsen T, Robakis TK, Bergink V, Mahjani B.
> **Genetic architecture of postpartum psychosis: from common to rare genetic variation.**
> *Molecular Psychiatry* (2026). https://doi.org/10.1038/s41380-026-03637-w

## Overview

This repository contains the code used to investigate the genetic architecture of postpartum psychosis (PP) through complementary analyses of common and rare genetic variation. Using family data from Swedish national registers and whole-genome sequencing (WGS) data from the *All of Us* Research Program, we:

- Estimated family-based heritability at **55%** (95% CI: 42–68%) using Falconer's Liability Threshold Model on sibling and cousin pairs.
- Estimated WGS-based heritability at **46%** (SE: 22%) using PCGC regression on common variants (MAF > 1%).
- Identified **HMGCR** as a high-confidence risk gene (FDR < 0.05) using the TADA model on 461 cases and 4,610 ancestry-matched controls.
- Replicated rare-variant burden associations for HMGCR and DNMT1 in the Mount Sinai BioMe Biobank (n = 58,990).
- Quantified genetic overlap between PP and bipolar disorder, schizophrenia, and 29 autoimmune diseases using the `propTrueNull` framework.

## Repository structure

```
pp_WGS/
├── heritability/          # WGS-based (PCGC) heritability
├── qc/                    # WGS sample- and variant-level QC (All of Us v7 + v8)
├── matching/              # PC-based case–control matching (optmatch)
├── annotation/            # VEP, LOFTEE, MPC-based variant classification
├── tada/                  # TADA gene-level association (Bayes factors, q-values)
├── replication/           # BioMe rare-variant replication (SKAT-O)
├── overlap/               # propTrueNull cross-disorder enrichment
└── figures/               # Manhattan plots, forest plots, histograms
```

## Requirements

- **R** (≥ 4.2): `data.table`, `optmatch`, `SKAT`, `limma`, `ggplot2`
- **Python** (≥ 3.9): `hail`, `pandas`, `numpy`, `scipy`
- **Tools**: VEP, LOFTEE, PLINK, PCGC
- Access to the *All of Us* Researcher Workbench (controlled tier)

## Data availability

- **All of Us** controlled tier data is available to registered institutions. Cohort selection details are in the workbench *"Detecting the prevalence of rare gene mutations."*
- **Swedish national register** data is restricted by law and available only via application to the Swedish Ethical Review Authority.
- **BioMe Biobank** data access is governed by Mount Sinai.


## Contact

For questions, please contact the corresponding author Behrang Mahjani (behrang.mahjani@mssm.edu) or open an issue.