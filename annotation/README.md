# Functional annotation and variant classification

Annotates ultra-rare variants with Ensembl VEP + LOFTEE and classifies them into functional categories (PTV, MisB, MisA, Mis0, SYN) used by the burden and TADA analyses.

## Pipeline

| Step | Script | Description |
|------|--------|-------------|
| 0 | 00_extract_ultrarare.py | Ultra-rare variants: AC ≤ 5 in cohort + gnomAD non-neuro |
| 1 | 01_prepare_vep_input.py | Export variants to VCF-style VEP input |
| 2 | 02_run_vep.sh | VEP (offline cache) + LOFTEE + Conservation plugins |
| 3 | 03_classify_variants.py | Classify into PTV / MisB / MisA / Mis0 / SYN |

> Prerequisite: cohort_matched.mt from matching/. VEP and LOFTEE resources (cache, FASTA, human_ancestor, GERP, conservation DB) must be installed on the host running step 02 — typically an HPC node, not the AoU workbench.

## Variant classification (paper Methods)

| Category | Definition |
|---|---|
| PTV | Protein-truncating (stop-gained, frameshift, splice donor/acceptor, transcript ablation); LOFTEE high-confidence |
| MisB | Missense, MPC ≥ 2 |
| MisA | Missense, 1 ≤ MPC < 2 |
| Mis0 | Missense, MPC < 1 |
| SYN | Synonymous (calibration set) |

Ultra-rare variants are defined as allele count ≤ 5 within the cohort and ≤ 5 (or absent) in the gnomAD v3.1 non-psychiatric subset.

## Software

- Ensembl VEP: https://www.ensembl.org/vep
- LOFTEE: https://github.com/konradjk/loftee
- MPC scores: Samocha et al. 2017
- Hail (≥ 0.2.130)

## References

- McLaren W et al. Genome Biol 2016 (VEP).
- Karczewski KJ et al. Nature 2020 (LOFTEE / gnomAD).
- Samocha KE et al. bioRxiv 2017 (MPC).
