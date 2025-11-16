# QLS600 Mendelian randomization tutorial

This repository contains R scripts and example data used in the QLS600 lecture on Mendelian randomization (MR).  
Related files are available for download at [OneDrive](https://mcgill-my.sharepoint.com/:f:/g/personal/satoshi_yoshiji_mcgill_ca/EqW78AwNRZJHialGi8OlcywBo84nzJwVYGNaMGGtAynDYg?e=CAqf8F) ..

The material walks through:

- A two-sample MR analysis of body mass index (BMI) on coronary artery disease (CAD)
- A protein-target MR and colocalization analysis of PCSK9 on CAD, using cis-pQTLs from ARIC

---

## Contents

- `BMI_noMHC_clumped.tsv`  
  Clumped genome-wide significant SNPs for BMI (exposure) outside the MHC region.  
  This file is formatted in the style returned by `TwoSampleMR::clump_data()` and is used as the exposure instrument set in `mr_bmi_to_cad.R`.

- `mr_bmi_to_cad.R`  
  Script demonstrating a basic two-sample MR pipeline:
  - Load BMI instruments (`BMI_noMHC_clumped.tsv`)
  - Load CAD GWAS summary statistics (`CAD/GCST90132314_buildGRCh38.formatted.tsv.gz`)
  - Harmonize exposure and outcome
  - Run standard MR methods (e.g. IVW, MR-Egger)
  - Run basic sensitivity analyses (heterogeneity, pleiotropy, Steiger directionality)
  - Save results and an MR scatter plot to the `output/` directory

- `mr_pcsk9_to_cad.R`  
  Script illustrating MR for a drug-like protein target and colocalization:
  - Exposure data: cis-pQTLs for PCSK9 (`PCSK9.5231_79.tsv`)
  - Outcome data: CAD GWAS (`CAD/GCST90132314_buildGRCh38.formatted.tsv.gz`)
  - Two-sample MR (PCSK9 → CAD) using `TwoSampleMR`
  - Sensitivity analyses (heterogeneity, pleiotropy, Steiger test)
  - Colocalization at the PCSK9 locus using `coloc::coloc.abf`, comparing:
    - Quantitative trait: circulating PCSK9 protein
    - Binary trait: CAD
  - Saves MR outputs, scatter plot, and colocalization summary into `output/`

**Note on PCSK9 data:**  
The PCSK9 cis-pQTL summary statistics (`PCSK9.5231_79.tsv`) used here are derived from the ARIC proteomics GWAS reported in *Nature Genetics* (PMID: **35501419**).

---

## Dependencies

Both scripts are written for R and rely on the following packages:

- [`TwoSampleMR`](https://mrcieu.github.io/TwoSampleMR/)
- [`tidyverse`](https://www.tidyverse.org/)
- [`data.table`](https://rdatatable.org/)
- [`coloc`](https://chr1swallace.github.io/coloc/) (for the PCSK9 colocalization section)

Install them, for example:

```r
install.packages(c("tidyverse", "data.table"))
# TwoSampleMR and coloc are usually installed from GitHub / Bioconductor:
# remotes::install_github("MRCIEU/TwoSampleMR")
# install.packages("coloc")
