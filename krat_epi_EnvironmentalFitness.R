# Final code for KratEpi Environmental Fitness Manuscript
# Lamka & Willoughby

# set working directory and read in dataset
setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/krat_epi/krat_epi/krat_epi_analysis/bwameth/Output") 
directory = getwd()
data = read.table("DATASET.csv", header=T, sep=",")

# below datasets taken from Harder et al. 2023
temp = read.table("temperature_mean.csv", header = T, sep=",") # temperature_max.csv
rain = read.table("PRISM_Rucker annual rain_forinport.csv", header = T, sep = ",")

# The goal is to find low inter-annual variation in methylation profiles, so will need to remove
# years with only one individual represented. We also need to remove siblings, as genetic relatedness
# may impact the variation in methylation profiles.

# remove individuals that are the only ones born in that year
table(data$year)
# years 1999 (indv 10), 2003, 2005 need to be removed
data = data[data$year != 1999 & data$year != 2003 & data$year != 2005,,drop=FALSE]

# remove siblings born in the same year to reduce confounding genetic relatedness
# from the pedigree we know that siblings include:
  # 3 & 4 are full sibs, 3 & 8 and 4 & 8 are half sibs, born 1997 -- two need to be removed
  # 35 & 36 are half sibs, born 1991 -- one needs to be removed
  # 7 & 29 are half sibs, born 1997 & 1993, so no worries here
  # 16 & 36 are half sibs, born 1992 & 1991, so no worries here
sibs_1997 <- c(3,4,8)
sample(sibs_1997, 1) # USE INDV 4
sibs_1991 <- c(35,36)
sample(sibs_1991, 1) # USE INDV 35
data = data[data$indvs != 3 & data$indvs != 8 & data$indvs != 36,,drop=FALSE]

# Calculate the coefficient of variation (CV = SD/mean) for each window
methylation_vars <- data %>%
  select(year, 35:14818) %>%  # Select year and methylation data
  pivot_longer(cols = -year, names_to = "window", values_to = "methylation") %>%
  group_by(year, window) %>%
  summarise(variation = sd(methylation, na.rm = TRUE) / mean(methylation, na.rm = TRUE), .groups = "drop")

# Find the lowest 0.05% of variation for each year
yearly_thresholds <- methylation_vars %>%
  group_by(year) %>%
  summarise(threshold = quantile(variation, probs = 0.005, na.rm = TRUE), .groups = "drop")

# Join back to methylation_vars to apply per-year thresholds
methylation_vars <- methylation_vars %>%
  left_join(yearly_thresholds, by = "year") %>%
  mutate(below_threshold = variation <= threshold)  # Flag windows below the 0.05% cutoff

# Filter for windows with the lowest variance
low_variance_windows <- methylation_vars %>%
  filter(below_threshold)

low_variance_windows # 888 total, 74 per year

range(low_variance_windows$threshold) # 0.0005303116 - 0.0833175819
mean(low_variance_windows$threshold) # 0.02060666
table(low_variance_windows$threshold,low_variance_windows$year)
table(data$year) # mean(table(data$year)) = 3.5
# note that threshold values are directly related to the number of indvs sampled in that year

# find windows that are signficant for more than one year
duplicates <- low_variance_windows %>%
  group_by(window) %>%
  filter(n() > 1) %>%
  ungroup()
# 854 unique windows, 31 windows in 2 years, 1 window in 4 years
# only 6 of the 32 duplicated windows were in successive years

# These windows with low inter-annual variation are considered the environmentally
# regulated windows in this system. Now, prep them for identifying their location
# along the banner-tailed kangaroo rat genome.

setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/krat_epi/krat_epi/krat_epi_analysis/bwameth/Output") 
filt <- read.table("filt_data_82.csv", header=T, sep=",") # grab raw dataset with the contig locations for each window

# Rename the windows and extract location information
v_columns = grep("^V", names(data), value = TRUE) 
short <- matrix(v_columns, nrow = nrow(filt), ncol = 1)
winds <- cbind(short, filt)
window_names_82 = cbind(winds$short, winds$contig, winds$start, winds$end)
colnames(window_names_82) <- c("short", "contig", "start", "end")

# Grab the name and location of the environmentally-regulated windows
envir_windows <- window_names_82[window_names_82[,1] %in% low_variance_windows$window, ]

nrow(envir_windows) # 854

# NOTE -- as of script creation, this was not updated and original .csv file was used -- not updated as of Oct 14, 2025
## write.csv(envir_windows, "signif_environmental_windows_nosibs_contig_start_end_june2025.csv", row.names = FALSE)
## "signif_environmental_windows_nosibs_contig_start_end_june2025.csv" == envir_windows

### DO THIS IN EASLEY

setwd("/home/gfl0003/kratepi/work/w_500_250")

# Load required packages
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(readr)

# Load significant methylation windows
sig_windows = read.table("signif_environmental_windows_nosibs_contig_start_end_june2025.csv", header=T, sep = ",")  
# Ensure it has columns: "contig", "start", "end"

# The windows are 500 bp with 250 bp overlap. We now need to add 1000 bp to each end of the window
# to include promoter regions. This is essentially because we have a contig-level genome assembly 
# and we want to make sure we include any important location information, as methylation in gene
# bodies and promoters have different functions.

for(x in 1:nrow(sig_windows)){
  sig_windows$start[x] <- sig_windows$start[x] - 1000
  sig_windows$end[x] <- sig_windows$end[x] + 1000
}

# Load GFF file
gff_file <- "GCF_019054845.1_ASM1905484v1_genomic.gff"
gff_data <- import(gff_file, format = "gff")

# Convert significant windows to GRanges
sig_gr <- GRanges(seqnames = sig_windows$contig,
                  ranges = IRanges(start = sig_windows$start, end = sig_windows$end))

# Find overlaps
overlaps <- findOverlaps(sig_gr, gff_data)

# Extract overlapping annotation details
overlapping_annotations <- data.frame(
  chr = seqnames(gff_data)[subjectHits(overlaps)],
  start = start(gff_data)[subjectHits(overlaps)],
  end = end(gff_data)[subjectHits(overlaps)],
  feature = gff_data$type[subjectHits(overlaps)],  # Feature type (gene, exon, etc.)
  gene_id = mcols(gff_data)$ID[subjectHits(overlaps)]  # Adjust based on your GFF annotation
)

# Check if 'Dbxref' exists and extract Entrez ID safely
if ("Dbxref" %in% colnames(mcols(gff_data))) {
  entrez_ids <- mcols(gff_data)$Dbxref[subjectHits(overlaps)]
  
  # Ensure consistency by handling cases where Dbxref is a list
  entrez_ids <- sapply(entrez_ids, function(x) {
    if (is.null(x) || length(x) == 0) {
      return(NA)  # Handle missing values
    } 
    if (any(grepl("GeneID:", x))) {
      return(sub(".*GeneID:([0-9]+).*", "\\1", paste(x, collapse = ";")))  # Extract numeric Entrez ID
    } 
    return(NA)  # If no GeneID found, return NA
  })
  
  # Add the Entrez ID column
  overlapping_annotations$entrez <- entrez_ids
} else {
  overlapping_annotations$entrez <- NA  # If 'Dbxref' doesn't exist, fill with NA
}

# Merge with significant windows
annotated_windows <- cbind(sig_windows[queryHits(overlaps), ], overlapping_annotations)

table(annotated_windows$feature) # all 854 windows aligned with annotated regions - great!

# Remove features that are just "regions", as those are not informative
annotated_windows <- annotated_windows[annotated_windows$feature != "region",,drop=FALSE] 

# Save results
write_csv(annotated_windows, "annotated_significant_environmental_windows_promoters_rmregion_june2025.csv")
# NOTE -- at the time of script creation, this .csv file was already on easley and used. This has not been
# updated as of Oct 14, 2025.

#### NOW MOVE BACK TO LOCAL COMPUTER

# Now, we want to grab only environmental windows and use those for downstream analyses. 
# Grab environmental windows from the dataset and intersect with the metadata

envir_windows_df = as.data.frame(envir_windows)
meta = data[,1:34]
methyl_windows = data[,35:ncol(data)]

# Ensure all values in envir_windows_df$short exist in methyl_windows column names
matching_columns <- colnames(methyl_windows) %in% envir_windows_df$short

# Subset methyl_windows by selecting only matching columns
methyl_subset <- methyl_windows[, matching_columns, drop = FALSE]
dim(methyl_subset)
# Now, we have the 854 environmentally regulated windows for all 42 individuals
# this will now be used for the Mr.PP / ridgeregression

### Prep abiotic factors for inclusion

# Prepare the environmental metadata for downstream analysis
temp = read.table("temperature_mean.csv", header = T, sep=",") # temperature_max.csv
rain = read.table("PRISM_Rucker annual rain_forinport.csv", header = T, sep = ",")

# Add rain$Mean to meta$rain based on matching year
meta$rain <- rain$ppt_inches[match(meta$year, rain$Date)] # SUM of rainfall

# Add temp$Mean to meta$temp based on matching year
meta$temp <- temp$Mean[match(meta$year, temp$Year)] # mean of monthly means

# added 10/14/25 -- add the march-march mean temp
meta$temp_March <- temp$Mean_March[match(meta$year, temp$Year)] # mean of monthly means

# now, grab population sizes
setwd("~/Library/CloudStorage/Box-Box/New Computer/Auburn/Data/KRats") 
krat = read.table("krat_pedinfo_years.csv", header=T, sep=",")

# Loop over each unique year in krat
for (yr in unique(krat$year)) {
  # Total population size for this year
  tot_size <- length(unique(krat$id[krat$year == yr]))
  
  # Assign to all rows in meta for that year
  meta$totpopsz[meta$year == yr] <- tot_size
  
  # # Now loop over subpops within this year -- need to change to numeric if want this
  # for (pop in unique(krat$subpop[krat$year == yr])) {
  #   subpop_size <- length(unique(krat$id[krat$year == yr & krat$subpop == pop]))
  #   
  #   meta$subpopsz[meta$year == yr & meta$pop == pop] <- subpop_size
  # }
}

remove(rain, temp, krat)
# NOTE that mean temp from march - march is curently added but this might want to be removed
# hard to find rainfall and pop size values for march-march, so might not be worth the new analyses

sig_methyl_data = cbind(meta, methyl_subset)
dim(sig_methyl_data) # 42 892 -- 38 meta + 854 environmental windows

########### REMOVE LATER #########
# See how response variables relate to environment
plot(data$year, data$off_survive)
plot(data$year, data$longevity)
points(data$year, data$off_survive, col = "red")
plot(temp$Year, temp$Mean, col = "green")

meta_cor <- cor(meta)
library(corrplot)
corrplot(meta_cor, method = "circle")
cor.test(meta$off_survive, meta$longevity) # cor = .57 ; signif
cor.test(meta$off_survive, meta$year) # cor = -.30 ; signif
cor.test(meta$longevity, meta$year) # cor = -.27 ; not signif
cor.test(meta$rain, meta$longevity) # cor = .18 ; not signif
cor.test(meta$temp, meta$longevity) # cor = -.27 ; not signif
cor.test(meta$rain, meta$off_survive) # cor = .33 ; signif
cor.test(meta$temp, meta$off_survive) # cor = -.33 ; signif

# Mr.PP
library(vegan)
# Create an empty data frame to store results
mrpp_results <- data.frame(Comparison = character(0),
                           A_value = numeric(0),
                           p_value = numeric(0),
                           stringsAsFactors = FALSE)

# Define MRPP comparisons
comparisons <- list(
  generation = sig_methyl_data$generation,
  lineage = sig_methyl_data$lineage,
  sex = sig_methyl_data$sex,
  genback = sig_methyl_data$genback,
  pop = sig_methyl_data$pop,
  month = sig_methyl_data$month,
  year = sig_methyl_data$year,
  age = sig_methyl_data$age,
  location = sig_methyl_data$location_numeric,
  disperser = sig_methyl_data$disp,
  temp = sig_methyl_data$temp,
  rain = sig_methyl_data$rain,
  totpopsz = sig_methyl_data$totpopsz
) 

# Loop through all MRPP tests and store results
for (comp in names(comparisons)) {
  mrpp_out <- vegan::mrpp(sig_methyl_data, grouping = comparisons[[comp]], 
                          permutations = 999, distance = "euclidean")
  
  # Store results in the data frame
  mrpp_results <- rbind(mrpp_results, 
                        data.frame(Comparison = comp, 
                                   A_value = mrpp_out$A, 
                                   p_value = mrpp_out$Pvalue))
}

# Pop size, rainfall, temp, year, lineage, subpop, month

# NOTE THIS HAS NOT BEEN UPDATED AS OF 10/15/25 -- NUMBERS NEED TO CHANGE SLIGHTLY IN TEXT -- PATTERNS ARE SAME
# write.csv(mrpp_results, "mrpp_results.csv", row.names = FALSE)

# NOTE these values are slightly off from what we reported before
# need to figure out why
# perhaps the envir windows arent aligned with the easley versions? need to check before submission!!!
# values are listed on kratepi_EnvironmentFitness_V1.docx
# code generated for the figs by Janna on lolipops.R (downloaded from Discord)

# NOTE - the only numbers I changed in the manuscript so far is the sample size per year. 
# one year had n = 10, mean = 3.75 in text
# I changed it to n=9, mean = 3.5
# this is likely due to the removed indviduals -- need to play with this to see which matches with downstream values


## RIDGE REGRESSIONS
library(ridge)
library(glmnet)
library(lme4)
library(caret)

# OFF_SURVIVE ~ MR.PP (pop size, rainfall, temp, year, lineage, subpop, month)
{
  combine_data_RRS <- data.frame(methyl_subset, sig_methyl_data[,c("lineage", "pop", "totpopsz",
                                                                   "month", "year", "rain", "temp",
                                                                   "off_survive")])
  
  # below suggested by chatGPT: (increase range of lambda and decrease repeats for faster identification)
  RRS <- caret::train(off_survive ~ .,
                      data = combine_data_RRS,
                      method = "glmnet",
                      tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(2, 3, length = 50)),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 5)) # 471
  
  # rr1 <- caret::train(off_survive ~ .,
  #                     data = combine_data_RRS,
  #                     method = "glmnet",
  #                     tuneGrid = expand.grid(alpha=0,lambda=seq(155,180,by=1)), 
  #                     trControl=trainControl(method="repeatedcv",n=10,repeats=100)) 
  
  best_RRS <- RRS$bestTune$lambda #the best lambda, need to  plug into next lines of code
  rr_RRS.out = linearRidge(off_survive ~ .,
                        data = combine_data_RRS,
                        lambda=best_RRS)
  summary(rr_RRS.out)
  
  # Extract coefficients summary
  ridge_summary_RRS <- summary(rr_RRS.out)
  
  # Extract p-values from the summary output
  pvals_RRS <- ridge_summary_RRS$summaries$summary1$coefficients[,5] # 5th column is Pr(>|t|)
  
  sigp_RRS <- pvals_RRS[pvals_RRS < 0.05]
  length(sigp_RRS) # 125

} # month, rain, totpopsiz

## For residual calcs -- 
# OFF_SURVIVE ~ METH
{
  off_survive <- sig_methyl_data[,8]
  combine_data_RRS_meth <- data.frame(methyl_subset, off_survive)
  
  # below suggested by chatGPT: (increase range of lambda and decrease repeats for faster identification)
  RRS_meth <- caret::train(off_survive ~ .,
                      data = combine_data_RRS_meth,
                      method = "glmnet",
                      tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(2, 3, length = 50)),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 5)) # 471
  
  # rr1 <- caret::train(off_survive ~ .,
  #                     data = combine_data_RRS,
  #                     method = "glmnet",
  #                     tuneGrid = expand.grid(alpha=0,lambda=seq(155,180,by=1)), 
  #                     trControl=trainControl(method="repeatedcv",n=10,repeats=100)) 
  
  best_RRS_meth <- RRS_meth$bestTune$lambda #the best lambda, need to  plug into next lines of code
  rr_RRS_meth.out = linearRidge(off_survive ~ .,
                           data = combine_data_RRS,
                           lambda=best_RRS_meth)
  summary(rr_RRS_meth.out)
  
  # Extract coefficients summary
  ridge_summary_RRS_meth <- summary(rr_RRS_meth.out)
  
  # Extract p-values from the summary output
  pvals_RRS_meth <- ridge_summary_RRS_meth$summaries$summary1$coefficients[,5] # 5th column is Pr(>|t|)
  
  sigp_RRS_meth <- pvals_RRS_meth[pvals_RRS_meth < 0.05]
  length(sigp_RRS_meth) # 125
  
} 

# OFF_SURVIVE ~ GENETICS (lineage, subpop)
{
  combine_data_RRS_gen <- data.frame(sig_methyl_data[,c("lineage", "pop", 
                                                                   "off_survive")])
  
  # below suggested by chatGPT: (increase range of lambda and decrease repeats for faster identification)
  RRS_gen <- caret::train(off_survive ~ .,
                      data = combine_data_RRS_gen,
                      method = "glmnet",
                      tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(0, 2, length = 50)),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 5)) # 471
  
  # rr1 <- caret::train(off_survive ~ .,
  #                     data = combine_data_RRS,
  #                     method = "glmnet",
  #                     tuneGrid = expand.grid(alpha=0,lambda=seq(155,180,by=1)), 
  #                     trControl=trainControl(method="repeatedcv",n=10,repeats=100)) 
  
  best_RRS_gen <- RRS_gen$bestTune$lambda #the best lambda, need to  plug into next lines of code
  rr_RRS_gen.out = linearRidge(off_survive ~ .,
                           data = combine_data_RRS_gen,
                           lambda=best_RRS_gen)
  summary(rr_RRS_gen.out)
  
  # Extract coefficients summary
  ridge_summary_RRS_gen <- summary(rr_RRS_gen.out)
  
  # Extract p-values from the summary output
  pvals_RRS_gen <- ridge_summary_RRS_gen$summaries$summary1$coefficients[,5] # 5th column is Pr(>|t|)
  
  sigp_RRS_gen <- pvals_RRS_gen[pvals_RRS_gen < 0.05]
  length(sigp_RRS_gen) # 125
  
}

# OFF_SURVIVE ~ ENVIR (pop size, rainfall, temp, year, month)
{
  combine_data_RRS_env <- data.frame(methyl_subset, sig_methyl_data[,c("totpopsz",
                                                                   "month", "year", "rain", "temp",
                                                                   "off_survive")])
  
  # below suggested by chatGPT: (increase range of lambda and decrease repeats for faster identification)
  RRS_env <- caret::train(off_survive ~ .,
                      data = combine_data_RRS_env,
                      method = "glmnet",
                      tuneGrid = expand.grid(alpha = 0, lambda = 10^seq(2, 3, length = 50)),
                      trControl = trainControl(method = "repeatedcv", number = 10, repeats = 5)) # 471
  
  # rr1 <- caret::train(off_survive ~ .,
  #                     data = combine_data_RRS,
  #                     method = "glmnet",
  #                     tuneGrid = expand.grid(alpha=0,lambda=seq(155,180,by=1)), 
  #                     trControl=trainControl(method="repeatedcv",n=10,repeats=100)) 
  
  best_RRS_env <- RRS_env$bestTune$lambda #the best lambda, need to  plug into next lines of code
  rr_RRS_env.out = linearRidge(off_survive ~ .,
                           data = combine_data_RRS_env,
                           lambda=best_RRS_env)
  summary(rr_RRS_env.out)
  
  # Extract coefficients summary
  ridge_summary_RRS_env <- summary(rr_RRS_env.out)
  
  # Extract p-values from the summary output
  pvals_RRS_env <- ridge_summary_RRS_env$summaries$summary1$coefficients[,5] # 5th column is Pr(>|t|)
  
  sigp_RRS_env <- pvals_RRS_env[pvals_RRS_env < 0.05]
  length(sigp_RRS_env) # 125
  
}










# code from ChatGPT
library(ggplot2)
library(dplyr)
library(tidyr)
library(lme4)
library(gridExtra)
library(tidyverse)
library(VennDiagram)

# Ridge Regression
library(ridge)
library(glmnet)
library(lme4)
library(caret)


  
  # # Calculate the overall mean and SD of variation across all windows and years
  # overall_mean <- mean(methylation_vars$variation, na.rm = TRUE)
  # overall_sd <- sd(methylation_vars$variation, na.rm = TRUE)
  # 

  # 
  # # Create a new dataframe for environmentally influenced windows -- want this from data, not methylation_variation. 
  # # PLUS, need to find which windows are consistently low for EVERY year
  # env_windows <- methylation_vars %>%
  #   filter(variation < threshold)
  # 
  # # Count number of low-variation windows per year
  # low_variation_counts <- env_windows %>%
  #   group_by(year) %>%
  #   summarise(num_windows = n(), .groups = "drop")
  # 
  # # Print the table
  # print(low_variation_counts)
  
  library(dplyr)
  

  
  # 0.1% has 270 (266 unique) -- 15 - 46 per year ?? uhh how?
  # ** 0.5% has 888 (856 unique) - 74 per year **
  # 1% has 1776 (1666 unqiue) -- 148 per year
  # 5% has 8881 (6607 unique) -- 740 per year
  # 10% has 17748 (10,130 unique) -- 1479 per year
  
  # # Compute mean and SD of variation for each year
  # yearly_stats <- methylation_vars %>%
  #   group_by(year) %>%
  #   summarise(yearly_mean = mean(variation, na.rm = TRUE),
  #             yearly_sd = sd(variation, na.rm = TRUE),
  #             .groups = "drop")
  # 
  # # Join back to methylation_vars to apply per-year thresholds
  # methylation_vars <- methylation_vars %>%
  #   left_join(yearly_stats, by = "year") %>%
  #   mutate(threshold = yearly_mean - 1 * yearly_sd)  # Adjust SD multiplier as needed
  
  # # Identify windows that are low in variation for **each specific year**
  # env_windows <- methylation_vars %>%
  #   filter(variation < threshold) %>%
  #   select(year, window, variation)  # Keep relevant columns
  # 
  # # Count number of low-variation windows per year
  # low_variation_counts <- env_windows %>%
  #   group_by(year) %>%
  #   summarise(num_windows = n(), .groups = "drop")
  # 
  # # Print the table
  # print(low_variation_counts)
  # 
  # sig_env_windows <- env_windows
  
  # Optionally, save the table as a CSV file
  # write.csv(sig_env_windows, "low_variation_windows_by_year_environment.csv", row.names = FALSE)
}

env = read.csv("low_variation_windows_by_year_environment.csv", header=T, sep=",")
