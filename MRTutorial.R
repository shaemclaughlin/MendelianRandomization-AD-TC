# Install packages
install.packages(c("devtools", "remotes", "knitr", "rmarkdown"))
remotes::install_github("MRCIEU/TwoSampleMR", force = TRUE)
library(TwoSampleMR)

# Set file paths to datasets
AD_file_path <- "/Users/shaemariemclaughlin/Desktop/MR_Data/AD.tsv"
TC_file_path <- "/Users/shaemariemclaughlin/Desktop/MR_Data/TC.tsv"

# Read exposure file into data.frame
TC_file <- read.csv(file = TC_file_path, sep = "\t")

# Read outcome file into data.frame
AD_file <- read.csv(file = AD_file_path, sep = "\t")

# Check if exposure object is data frame
is.data.frame(TC_file)

# Check if outcome object is data frame
is.data.frame(AD_file)

# Print colnames for exposure and outcome data frames
colnames(AD_file)
colnames(TC_file)

# Convert existing exposure data frame into correct format using format_data
exp_dat <- format_data(TC_file, type = "exposure", snps = NULL, header = TRUE, phenotype_col = "TRAIT", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "AF", effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "P", z_col = "Z", chr_col = "CHROM", pos_col = "POS", samplesize_col = "N")
exp_dat
any(is.na(exp_dat))
exp_dat <- na.omit(exp_dat)
any(is.na(exp_dat))
nrow(exp_dat)

n_splits <- 10
chunk_size <- ceiling(nrow(exp_dat)/n_splits)
clumped_exp_dat_list <- lapply(seq_len(n_splits), function(i){
  start_row <- (i - 1) * chunk_size + 1
  end_row <- min(i * chunk_size, nrow(exp_dat))
  chunk_dat <- exp_dat[start_row:end_row,]
  clump_data(chunk_dat)
})
clumped_exp_dat <- do.call(rbind, clumped_exp_dat_list)
  
nrow(clumped_exp_dat)
any(is.na(clumped_exp_dat))

# Convert existing outcome data frame into correct format using format_data
outcome_dat <- format_data(AD_file, type = "outcome", snps = NULL, header = TRUE, phenotype_col = "TRAIT", snp_col = "SNP", beta_col = "BETA", se_col = "SE", eaf_col = "AF", effect_allele_col = "ALT", other_allele_col = "REF", pval_col = "P", z_col = "Z", chr_col = "CHROM", pos_col = "POS", samplesize_col = "N")
outcome_dat
any(is.na(outcome_dat))
outcome_dat <- na.omit(outcome_dat)
any(is.na(outcome_dat))
nrow(outcome_dat)

# Test on a smaller datasets (first 100,000 rows of exposure and outcome data)
thousand_row_exp_dat <- exp_dat[1:200000,]
thousand_row_out_dat <- outcome_dat[1:200000,]

# Harmonize the exposure and outcome data
harmonized_dat <- harmonise_data(clumped_exp_dat, outcome_dat)
nrow(harmonized_dat)

# Check harmonized data for NAs
any(is.na(harmonized_dat))
nrow(harmonized_dat)
harmonized_dat <- na.omit(harmonized_dat)
nrow(harmonized_dat)

harmonized_dat$Beta_IV <- harmonized_dat$beta.outcome/harmonized_dat$beta.exposure
hist(harmonized_dat$Beta_IV, breaks=50, main="Beta_IV")
any(is.na(harmonized_dat$Beta_IV))
str(harmonized_dat)
summary(harmonized_dat$Beta_IV)

# Identify infinite values
inf_rows <- is.infinite(harmonized_dat$Beta_IV)
sum(inf_rows)
harmonized_dat[inf_rows,]

harmonized_dat <- harmonized_dat[!inf_rows,]
any(is.infinite(harmonized_dat$Beta_IV))
summary(harmonized_dat$Beta_IV)

# Perform MR analysis
res <- mr(harmonized_dat)
res

p1 <- mr_scatter_plot(res, harmonized_dat)
p1[[1]]
length(p1)

het_res <- mr_heterogeneity(harmonized_dat)
het_res

hp_res <- mr_pleiotropy_test(harmonized_dat)
hp_res

res_single <- mr_singlesnp(harmonized_dat)
res_single
str(res_single)

res_loo <- mr_leaveoneout(harmonized_dat)
res_loo


