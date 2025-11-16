########################################
# Step 1 MR: BMI (exposure) -> CAD (outcome)
########################################

# -------------------------
# 1. Load required libraries
# -------------------------
# TwoSampleMR: core functions for MR (read_exposure_data, harmonise_data, mr, etc.)
# tidyverse: for convenient data manipulation (dplyr, readr, etc.)
# data.table: we use fread() to read large GWAS files efficiently

library(TwoSampleMR)
library(tidyverse)
library(data.table)

# -------------------------
# 2. Set working directory and create output folders
# -------------------------
# This is the folder where all input and output files live.
wd <- "/Users/sy/Documents/Studies/0.YoshijiLab/11.lecture/QLS600/MR_tutorial/"
setwd(wd)

# Create an "output" folder and subfolders for each type of result.
output_dir <- file.path(wd, "output/")
dirs <- c("harmonized", "or", "hetero", "pleio", "steiger")
for (dir in dirs) {
  system(paste0("mkdir -p ", output_dir, dir, "/"))
}

# Name for this MR analysis (used as prefix for output files)
protname <- "BMI_to_CAD"

# -------------------------
# 3. Define paths to exposure (BMI) and outcome (CAD) files
# -------------------------
# BMI exposure data: clumped SNPs with effect sizes, SE, allele frequencies, etc.
exp_path <- file.path(wd, "BMI_noMHC_clumped.tsv")

# CAD outcome GWAS: summary statistics with rs_id, beta, SE, etc.
outcome_path <- file.path(
  wd,
  "CAD/GCST90132314_buildGRCh38.formatted.tsv.gz"
)

# -------------------------
# 4. Read outcome (CAD) GWAS using fread()
# -------------------------
# CAD file columns are:
# chromosome, id, base_pair_position, effect_allele, other_allele,
# p_value, effect_allele_frequency, beta, standard_error, case, sample_size, rs_id
#
# We read with fread() (fast I/O), convert to data.frame,
# and rename rs_id to SNP, because TwoSampleMR expects a column called "SNP".
outcome_GWAS <- fread(outcome_path, data.table = FALSE)
outcome_GWAS <- dplyr::rename(outcome_GWAS, SNP = rs_id)

# -------------------------
# 5. Read exposure (BMI) data using TwoSampleMR helper
# -------------------------
# This file already has the column names in the right format:
# SNP, beta.exposure, se.exposure, effect_allele.exposure,
# other_allele.exposure, eaf.exposure, pval.exposure, samplesize.exposure
exp_dat <- read_exposure_data(
  filename          = exp_path,
  sep               = "\t",
  snp_col           = "SNP",
  beta_col          = "beta.exposure",
  se_col            = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col  = "other_allele.exposure",
  eaf_col           = "eaf.exposure",
  pval_col          = "pval.exposure",
  samplesize_col    = "samplesize.exposure"
)

# -------------------------
# 6. Restrict CAD GWAS to BMI instruments
# -------------------------
# We keep only SNPs that are present in the BMI instrument set.
outcome_GWAS <- dplyr::filter(outcome_GWAS, SNP %in% exp_dat$SNP)

# -------------------------
# 7. Format CAD data for TwoSampleMR as outcome
# -------------------------
# format_data() converts a generic GWAS data.frame into the standardized format
# that TwoSampleMR expects for outcome data.
formatted_outcome <- format_data(
  outcome_GWAS,
  snps              = exp_dat$SNP,
  type              = "outcome",
  snp_col           = "SNP",
  beta_col          = "beta",
  se_col            = "standard_error",
  eaf_col           = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col  = "other_allele",
  pval_col          = "p_value",
  chr_col           = "chromosome",
  pos_col           = "base_pair_position",
  samplesize_col    = "sample_size",
  ncase_col         = "case",
  id = "CAD"
)

# -------------------------
# 8. Harmonise BMI exposure and CAD outcome
# -------------------------
# harmonise_data() aligns effect alleles between exposure and outcome.
# We then filter on allele frequency to avoid very rare variants.
exp_dat_outcome <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat  = formatted_outcome
)

exp_dat_outcome <- dplyr::filter(
  exp_dat_outcome,
  (eaf.exposure > 0.001 & eaf.exposure < 0.999) &
    (eaf.outcome  > 0.001 & eaf.outcome  < 0.999)
)

# Save harmonised dataset for record and possible plotting later.
exp_data_outcome_name <- paste0(output_dir, "harmonized/", protname, ".harmonized.txt")
write.table(exp_dat_outcome, file = exp_data_outcome_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 9. Run main MR analyses and output odds ratios
# -------------------------
# mr() runs several MR methods (e.g. IVW, MR-Egger).
# generate_odds_ratios() converts log-odds estimates to odds ratios
# for easier interpretation when the outcome is binary (CAD here).
mr_results <- mr(exp_dat_outcome)
OR <- generate_odds_ratios(mr_results)

OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write.table(OR, file = OR_name,
            sep = "\t", quote = FALSE, row.names = FALSE)


# -------------------------
# 10. Create MR scatter plot (BMI vs CAD)
# -------------------------
# mr_scatter_plot() creates a scatter plot of SNP-exposure vs SNP-outcome effects
# with the fitted MR line. It returns a list of ggplot objects (one per method/outcome).
# For this simple example, we use the first plot and save it as a PNG.
scatter_list <- mr_scatter_plot(mr_results, exp_dat_outcome)

# Take the first scatter plot in the list
scatter_plot <- scatter_list[[1]]
print(scatter_plot)

# Save the scatter plot to the "plots" folder
scatter_file <- paste0(output_dir, "plots/", protname, "_scatter.png")
ggplot2::ggsave(filename = scatter_file,
                plot     = scatter_plot,
                width    = 6,
                height   = 6,
                dpi      = 300)

# -------------------------
# 11. MR-Egger intercept (horizontal pleiotropy) test
# -------------------------
# mr_pleiotropy_test() checks for directional pleiotropy using the MR-Egger intercept.
pleio_res <- mr_pleiotropy_test(exp_dat_outcome)

pleio_name <- paste0(output_dir, "pleio/", protname, ".pleio.txt")
write.table(pleio_res, file = pleio_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 12. Heterogeneity test (Cochran's Q and I^2)
# -------------------------
# mr_heterogeneity() computes heterogeneity statistics across SNPs.
# We also calculate I^2 as a simple measure of heterogeneity.
hetero_res <- mr_heterogeneity(exp_dat_outcome)
hetero_res$isquared <- abs(100 * (hetero_res$Q - hetero_res$Q_df) / hetero_res$Q)

hetero_name <- paste0(output_dir, "hetero/", protname, ".hetero.txt")
write.table(hetero_res, file = hetero_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 13. Steiger directionality test
# -------------------------
# directionality_test() checks whether the instruments explain more variance
# in the exposure than in the outcome (to support the direction BMI -> CAD).
# Here, samplesize.exposure comes from the BMI file and
# samplesize.outcome from the CAD file via format_data().
steiger <- directionality_test(exp_dat_outcome)

steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write.table(steiger, file = steiger_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

