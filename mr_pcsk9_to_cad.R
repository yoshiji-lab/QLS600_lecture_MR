########################################
# Step 2 MR: PCSK9 (exposure) -> CAD (outcome)
########################################

# -------------------------
# 1. Load required libraries
# -------------------------
# TwoSampleMR: core functions for MR (read_exposure_data, harmonise_data, mr, etc.)
# tidyverse: for convenient data manipulation and plotting (dplyr, ggplot2, readr, etc.)
# data.table: we use fread() to read large GWAS files efficiently

library(TwoSampleMR)
library(tidyverse)
library(data.table)

# -------------------------
# 2. Set working directory and create output folders
# -------------------------
# This is the folder where all input and output files live.
# Set your own wd!
wd <- "???"
setwd(wd)

# Create an "output" folder and subfolders for each type of result.
output_dir <- file.path(wd, "output/")
dirs <- c("harmonized", "or", "hetero", "pleio", "steiger",
          "plots")
for (dir in dirs) {
  system(paste0("mkdir -p ", output_dir, dir, "/"))
}

# Name for this MR analysis (used as prefix for output files)
protname <- "PCSK9_to_CAD"

# -------------------------
# 3. Define paths to exposure (PCSK9) and outcome (CAD) files
# -------------------------
# PCSK9 exposure data: cis pQTLs for PCSK9.5231_79
# Columns include:
# protein, seqid, rsid, Effect_allele, Other_allele,
# Effect_allele_freq, BETA, SE, P, N, ...
exp_path <- file.path(wd, "PCSK9.5231_79.tsv")

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
# 5. Read exposure (PCSK9) data using TwoSampleMR helper
# -------------------------
# PCSK9 file columns:
# rsid, Effect_allele, Other_allele, Effect_allele_freq, BETA, SE, P, N, ...
#
# We map these to the standard exposure fields.

#Look into the data to assign these columns!

exp_dat <- read_exposure_data(
  filename          = exp_path,
  sep               = "\t",
  snp_col           = "rsid",
  beta_col          = "???",
  se_col            = "???",
  effect_allele_col = "???",
  other_allele_col  = "Other_allele",
  eaf_col           = "???",
  pval_col          = "???",
  samplesize_col    = "???"
)

# -------------------------
# 6. Restrict CAD GWAS to PCSK9 instruments
# -------------------------
# We keep only SNPs that are present in the PCSK9 cis pQTL set.
outcome_GWAS <- dplyr::filter(outcome_GWAS, SNP %in% exp_dat$SNP)

# -------------------------
# 7. Format CAD data for TwoSampleMR as outcome
# -------------------------
# format_data() converts a generic GWAS data.frame into the standardized format
# that TwoSampleMR expects for outcome data.

#Look into the data to assign these columns!

formatted_outcome <- format_data(
  outcome_GWAS,
  snps              = exp_dat$SNP,
  type              = "outcome",
  snp_col           = "???",
  beta_col          = "???",
  se_col            = "???",
  eaf_col           = "???",
  effect_allele_col = "???",
  other_allele_col  = "???",
  pval_col          = "???",
  chr_col           = "???",
  pos_col           = "???",
  samplesize_col    = "???",
  ncase_col         = "???"
)

# -------------------------
# 8. Harmonise PCSK9 exposure and CAD outcome
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
# mr() runs several MR methods (for example IVW, MR-Egger).
# generate_odds_ratios() converts log-odds estimates to odds ratios
# for easier interpretation when the outcome is binary (CAD here).
mr_results <- mr(exp_dat_outcome)
OR <- generate_odds_ratios(???)

OR_name <- paste0(output_dir, "or/", protname, ".or.txt")
write.table(OR, file = OR_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 10. Create MR scatter plot (PCSK9 vs CAD)
# -------------------------
# mr_scatter_plot() creates a scatter plot of SNP-exposure vs SNP-outcome effects
# with the fitted MR line. It returns a list of ggplot objects.
scatter_list <- mr_scatter_plot(???, exp_dat_outcome)

# Take the first scatter plot in the list
# How do we specify the first one?
scatter_plot <- scatter_list[[???]]

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
# in the exposure than in the outcome. Here that means PCSK9 -> CAD.
steiger <- directionality_test(exp_dat_outcome)

steiger_name <- paste0(output_dir, "steiger/", protname, ".steiger.txt")
write.table(steiger, file = steiger_name,
            sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------
# 14. MR-PRESSO (outlier and distortion test)
# -------------------------
# MR-PRESSO detects horizontal pleiotropic outliers and tests whether
# omitted

#Skipping 14

# -------------------------
# 15. Final message
# -------------------------
print(paste0(protname, " MR analysis completed."))

# -------------------------------------------------------------------------------------

########################################
# 16. Colocalization: PCSK9 protein vs CAD
########################################

# -------------------------
# 16.1 Load coloc library
# -------------------------
# coloc: functions to perform colocalization analysis (e.g. coloc.abf).
library(coloc)

# -------------------------
# 16.2 Load PCSK9 protein GWAS (quantitative trait)
# -------------------------
# We treat PCSK9.5231_79.tsv as the protein GWAS for PCSK9.
# It contains cis-pQTLs with effect sizes, SE, p-values, and allele frequencies.
pcsk9_path <- file.path(wd, "PCSK9.5231_79.tsv")
pcsk9_gwas <- data.table::fread(pcsk9_path, data.table = FALSE)

# Rename columns to generic names expected by coloc:
#   snp      : variant ID (rsID)
#   beta     : effect size (BETA)
#   SE       : standard error
#   MAF      : minor/effect allele frequency
#   pvalues  : p-value
pcsk9_gwas <- dplyr::rename(
  pcsk9_gwas,
  snp     = rsid,
  beta    = ???,
  SE      = ???,
  MAF     = ???,
  pvalues = ???
)

# Compute variance of beta (varbeta = SE^2)
pcsk9_gwas$varbeta <- pcsk9_gwas$SE^2

# Clean up edge cases for MAF and p-values to avoid numerical issues
pcsk9_gwas$MAF[pcsk9_gwas$MAF == 0]    <- 1e-4
pcsk9_gwas$MAF[pcsk9_gwas$MAF == 1]    <- 0.9999
pcsk9_gwas$pvalues[pcsk9_gwas$pvalues == 0] <- 1e-300
pcsk9_gwas$pvalues[pcsk9_gwas$pvalues == 1] <- 0.9999

# Use a single sample size for the protein GWAS.
# We take the first unique N from the file (all rows should have the same N).
N_pcsk9 <- unique(pcsk9_gwas$N)[1]

# -------------------------
# 16.3 Load full CAD GWAS for colocalization
# -------------------------
# We reload the full CAD summary statistics (not restricted to MR instruments)
# because colocalization should use all overlapping SNPs in the locus.
cad_full <- data.table::fread(outcome_path, data.table = FALSE)

# Rename columns to generic names:
#   snp      : variant ID (rsID)
#   beta     : log-odds ratio for CAD
#   SE       : standard error
#   MAF      : effect allele frequency
#   pvalues  : p-value
cad_df <- dplyr::rename(
  cad_full,
  snp     = rs_id,
  beta    = ???,
  SE      = ???,
  MAF     = ???,
  pvalues = ???
)

# Compute variance of beta (varbeta = SE^2)
cad_df$varbeta <- cad_df$SE^2

# Clean up edge cases for MAF and p-values
cad_df$MAF[cad_df$MAF == 0]    <- 1e-4
cad_df$MAF[cad_df$MAF == 1]    <- 0.9999
cad_df$pvalues[cad_df$pvalues == 0] <- 1e-300
cad_df$pvalues[cad_df$pvalues == 1] <- 0.9999

# Derive total sample size (N) and case proportion (s) for CAD.
# Here we use the maximum values across rows as a simple approximation.
N_cad     <- max(cad_full$sample_size, na.rm = TRUE)
cases_cad <- max(cad_full$case,         na.rm = TRUE)
s_cad     <- cases_cad / N_cad

# -------------------------
# 16.4 Find overlapping SNPs between PCSK9 and CAD
# -------------------------
# We restrict colocalization to SNPs that are present in both datasets.
common <- dplyr::inner_join(
  pcsk9_gwas,
  cad_df,
  by      = "snp",
  suffix  = c(".pcsk9", ".cad")
)

# Check how many common SNPs we have and stop if there are none.
n_common <- nrow(common)
message("Number of common SNPs for PCSK9 vs CAD: ", n_common)
stopifnot(n_common > 0)

# -------------------------
# 16.5 Save common SNPs and prepare coloc inputs
# -------------------------
# All colocalization-related files go into output/coloc/
coloc_dir <- file.path(output_dir, "coloc")
system(paste0("mkdir -p ", coloc_dir))

# Save common SNPs for inspection
common_file <- file.path(coloc_dir, paste0(protname, "_PCSK9_CAD_common.tsv"))
write.table(
  common,
  file      = common_file,
  sep       = "\t",
  quote     = FALSE,
  row.names = FALSE
)

# Dataset 1: PCSK9 protein GWAS (quantitative trait)
sumstats1_list <- list(
  type    = "quant",
  N       = N_pcsk9,
  MAF     = common$MAF.pcsk9,
  beta    = common$beta.pcsk9,
  varbeta = common$varbeta.pcsk9,
  pvalues = common$pvalues.pcsk9
)

# Dataset 2: CAD GWAS (case-control)
sumstats2_list <- list(
  type    = "cc",
  N       = N_cad,
  s       = s_cad,
  MAF     = common$MAF.cad,
  beta    = common$beta.cad,
  varbeta = common$varbeta.cad,
  pvalues = common$pvalues.cad
)

# -------------------------
# 16.6 Run colocalization analysis
# -------------------------
# coloc.abf() computes posterior probabilities that:
#  H0: neither trait has a causal variant
#  H1: only trait 1 has a causal variant
#  H2: only trait 2 has a causal variant
#  H3: both have causal variants, different
#  H4: both share the same causal variant (colocalization)
coloc_res <- coloc::coloc.abf(
  dataset1 = sumstats1_list,
  dataset2 = sumstats2_list
)

print(coloc_res)
# -------------------------
# 16.7 Save colocalization summary
# -------------------------
# We extract the summary table and save it in output/coloc/
coloc_file <- file.path(coloc_dir, paste0(protname, "_coloc_summary.tsv"))

coloc_summary_df <- as.data.frame(coloc_res$summary)
coloc_summary_df <- tibble::rownames_to_column(coloc_summary_df, var = "hypothesis")

readr::write_tsv(coloc_summary_df, file = coloc_file)

message("Colocalization analysis (PCSK9 vs CAD) completed and saved to: ", coloc_file)

