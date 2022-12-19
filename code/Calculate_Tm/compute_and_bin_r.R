#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 3){stop("Usage: <tvec> <prefix for test panel genotypes> <outpath>") }

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
}))

tvec_file=args[1] # Test vec
test_file_prefix=args[2] # Test panel prefix
outpath = args[3]

# Read in test vec
Tvec <- fread(tvec_file)
n <- nrow(Tvec)

# Compute t(X)T using plink
outfile_XT <- paste0(outpath, "xt")
cmd_XT <- paste("sh code/Calculate_Tm/compute_XT.sh", test_file_prefix, tvec_file, outfile_XT, sep = " ")
system(cmd_XT)

# Adjust Betas to account for variance in x

# Read in betas and genotype counts
beta_plink <- fread(paste0(outpath, "xt.Tvec.glm.linear"))
count_plink <- fread(paste0(outpath, "xt.gcount"))

# Calculate length of mean centered genotypes from counts
nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
mean_gc <- counts / nOBS
length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

# Fix betas
betas_plink_norm <- beta_plink$BETA * length_mc_genos

#  Re-write .linear file with correct betas
beta_plink$BETA <- betas_plink_norm
beta_reformat <- beta_plink %>% dplyr::select(ID, A1, BETA)
beta_reformat$r <- beta_reformat$BETA / (n -1)

# Save r
fwrite(beta_reformat, paste0(outpath, "r.txt"),row.names=F,quote=F,sep="\t", col.names = T)

## Bin r based on deciles
bin_r <- beta_reformat %>% select(r) %>% mutate(r_bin = ntile(r, n=10))%>% group_by(r_bin) %>% summarise(avg = mean(r), SD = sd(r))

# Save bins
fwrite(bin_r, paste0(outpath, "r_bins.txt"),row.names=F,quote=F,sep="\t", col.names = T)
