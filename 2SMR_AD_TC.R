#Packages 
library(tidyverse)
library(TwoSampleMR)
library(rmarkdown)
library(RadialMR)

#Files 
AD_file_path <- "/Users/shaemariemclaughlin/Desktop/MR_Data/AD.tsv" 
TC_file_path <- "/Users/shaemariemclaughlin/Desktop/MR_Data/TC.tsv"

#The exposure in this analysis is total cholesterol 
#Formatting exposure data
exposure_data <- read.csv(file = TC_file_path, sep = "\t")
exposure_formatted <- format_data(exposure_data, 
                                  type ="exposure", 
                                  snps = NULL, 
                                  header = TRUE, 
                                  phenotype_col = "TRAIT", 
                                  snp_col = "SNP",
                                  beta_col = "BETA", 
                                  se_col = "SE", 
                                  eaf_col = "AF", 
                                  effect_allele_col = "ALT", 
                                  other_allele_col = "REF", 
                                  pval_col = "P", 
                                  z_col = "Z", 
                                  chr_col = "CHROM", 
                                  pos_col = "POS", 
                                  samplesize_col = "N")

#Clumping exposure data 
exposure_clumped <- exposure_formatted %>%
filter(pval.exposure < 1e-6) %>% clump_data(., )

#The outcome in this analysis is Alzheimer's Disease 
#Formatting outcome data
outcome_data <- read.csv(file = AD_file_path, sep = "\t")
outcome_formatted <- format_data(outcome_data, 
                                 type ="outcome", 
                                 snps = NULL, 
                                 header = TRUE, 
                                 phenotype_col = "TRAIT", 
                                 snp_col = "SNP", 
                                 beta_col = "BETA", 
                                 se_col = "SE", 
                                 eaf_col = "AF", 
                                 effect_allele_col = "ALT",
                                 other_allele_col = "REF", 
                                 pval_col = "P", 
                                 z_col = "Z", 
                                 chr_col = "CHROM", 
                                 pos_col = "POS", 
                                 samplesize_col = "N")

#Extracting outcome data 
outcome_clumped <- filter(outcome_formatted, SNP %in% exposure_clumped$SNP)

#Harmonization 
harmonized_data <- harmonise_data(exposure_dat = exposure_clumped, 
                                  outcome_dat = outcome_clumped)

#Performing Mendelian randomization 
results <- mr(harmonized_data,method_list = c("mr_egger_regression", 
                                              "mr_weighted_median",
                                              "mr_ivw_fe", 
                                              "mr_weighted_mode")) 
knitr::kable(results[,3:9], col.names = c("Outcome", 
                                          "Exposure", 
                                          "Method", 
                                          "NSNP", 
                                          "B", 
                                          "SE", 
                                          "P-Value"))

#Getting heterogeneity statistics 
heterogeneity_results <- mr_heterogeneity(harmonized_data) 
knitr::kable(heterogeneity_results[,3:8], col.names = c("Outcome", 
                                                        "Exposure", 
                                                        "Method", 
                                                        "Q", 
                                                        "Q (DF)", 
                                                        "Q P-Value"))

#Performing horizontal pleiotropy test 
hp_results <- mr_pleiotropy_test(harmonized_data) 
knitr::kable(hp_results[,3:7], col.names = c("Outcome", 
                                             "Exposure", 
                                             "Egger Intercept", 
                                             "SE", 
                                             "P-Value"))

#Performing leave one out Mendelian randomization 
loo_results <- mr_leaveoneout(harmonized_data)
knitr::kable(loo_results[,c(1, 2, 5, 6, 7, 8, 9)], col.names = c("Exposure", 
                                                                 "Outcome", 
                                                                 "Sample Size", 
                                                                 "SNP", 
                                                                 "B", 
                                                                 "SE", 
                                                                 "P-Value"))

#Performing single snp Mendelian randomization 
singlesnp_results <- mr_singlesnp(harmonized_data)
knitr::kable(singlesnp_results[,c(1, 2, 5, 6, 7, 8, 9)], col.names = c("Exposure", 
                                                                       "Outcome", 
                                                                       "Sample Size", 
                                                                       "SNP", 
                                                                       "B", 
                                                                       "SE", 
                                                                       "P-Value"))

#Setting plot theme 
theme_set(theme_bw())

#Scatter plot 
scatter <- mr_scatter_plot(results, harmonized_data)
scatter[[1]]

#Forest plot 
forest <- mr_forest_plot(singlesnp_results) 
forest[[1]]

#Funnel plot 
funnel <- mr_funnel_plot(singlesnp_results) 
funnel[[1]]

#Leave-one-out plot 
loo <- mr_leaveoneout_plot(loo_results) 
loo[[1]]

#Radial formatting
radial_data <- harmonized_data %>% filter(mr_keep == TRUE) %>% dat_to_RadialMR()
str(radial_data)
radial_mr_data <- radial_data$TC.AD

#Radial IVW
bonff = 0.05/nrow(radial_mr_data) #Bonferroni correction to account for multiple testing
radial_ivw_results <- ivw_radial(radial_mr_data, alpha = bonff)
radial_ivw_results_data <- radial_ivw_results$data

#Getting radial IVW outlier snps
radial_ivw_outlier_snps <- radial_ivw_results$outliers

#outlier snp table
knitr::kable(radial_ivw_outlier_snps[], col.names = c("Outlier SNPs", "Q Statistic", "P-Value"))

#Radial Egger
radial_egger_results <- egger_radial(radial_mr_data, alpha = bonff)
radial_egger_results
radial_egger_results_data <- radial_egger_results$data


#Getting radial Egger outlier snps
radial_egger_outlier_snps <- radial_egger_results$outliers

#outlier snp table
knitr::kable(radial_egger_outlier_snps[], col.names = c("Outlier SNPs", "Q Statistic", "P-Value"))

#Radial plots
ivw_radial_plot <- plot_radial(radial_ivw_results, radial_scale = FALSE, show_outliers = FALSE)

egger_radial_plot <- plot_radial(radial_egger_results, radial_scale = FALSE, show_outliers = FALSE)

cowplot::plot_grid(
  ivw_radial_plot + coord_fixed(ratio=0.25) + theme(legend.position = 'bottom'), 
  egger_radial_plot + theme(legend.position = 'bottom'), 
  align = 'h'
)

#joining outlier snps
combined_outliers <- bind_rows(radial_egger_outlier_snps, radial_ivw_outlier_snps) %>% distinct(SNP, .keep_all = TRUE)

#removing outlier snps from harmonized data
harmonized_data_no_outliers <- harmonized_data[!harmonized_data$SNP %in% c("rs10468017", 
                                                                           "rs1883025", 
                                                                           "rs515135", 
                                                                           "rs6504872", 
                                                                           "rs6544713", 
                                                                           "rs7412", 
                                                                           "rs75687619", 
                                                                           "rs8103315"),]

#Performing Mendelian randomization with outliers removed
results <- mr(harmonized_data_no_outliers,method_list = c("mr_egger_regression", 
                                              "mr_weighted_median",
                                              "mr_ivw_fe", 
                                              "mr_weighted_mode")) 
knitr::kable(results[,3:9], col.names = c("Outcome", 
                                          "Exposure", 
                                          "Method", 
                                          "NSNP", 
                                          "B", 
                                          "SE", 
                                          "P-Value"))

#Getting heterogeneity statistics 
heterogeneity_results <- mr_heterogeneity(harmonized_data_no_outliers) 
knitr::kable(heterogeneity_results[,3:8], col.names = c("Outcome", 
                                                        "Exposure", 
                                                        "Method", 
                                                        "Q", 
                                                        "Q (DF)", 
                                                        "Q P-Value"))

#Performing horizontal pleiotropy test 
hp_results <- mr_pleiotropy_test(harmonized_data_no_outliers) 
knitr::kable(hp_results[,3:7], col.names = c("Outcome", 
                                             "Exposure", 
                                             "Egger Intercept", 
                                             "SE", 
                                             "P-Value"))

#Performing leave one out Mendelian randomization 
loo_results <- mr_leaveoneout(harmonized_data_no_outliers)
knitr::kable(loo_results[,c(1, 2, 5, 6, 7, 8, 9)], col.names = c("Exposure", 
                                                                 "Outcome", 
                                                                 "Sample Size", 
                                                                 "SNP", 
                                                                 "B", 
                                                                 "SE", 
                                                                 "P-Value"))

#Performing single snp Mendelian randomization 
singlesnp_results <- mr_singlesnp(harmonized_data_no_outliers)
knitr::kable(singlesnp_results[,c(1, 2, 5, 6, 7, 8, 9)], col.names = c("Exposure", 
                                                                       "Outcome", 
                                                                       "Sample Size", 
                                                                       "SNP", 
                                                                       "B", 
                                                                       "SE", 
                                                                       "P-Value"))

#Setting plot theme 
theme_set(theme_bw())

#Scatter plot 
scatter <- mr_scatter_plot(results, harmonized_data_no_outliers)
scatter[[1]]

#Forest plot 
forest <- mr_forest_plot(singlesnp_results) 
forest[[1]]

#Funnel plot 
funnel <- mr_funnel_plot(singlesnp_results) 
funnel[[1]]

#Leave-one-out plot 
loo <- mr_leaveoneout_plot(loo_results) 
loo[[1]]

