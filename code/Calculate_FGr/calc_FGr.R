# This script computes FGr by doing GX'T

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_FGr.R <test_prefix> <gwas_prefix> <tvec_file> <out_prefix>")}

suppressWarnings(suppressMessages({
  library(pgenlibr)
  library(data.table)
  library(dplyr)
  library(pracma)
}))

test_prefix = args[1] # Prefix to GWAS plink files
#test_prefix = "output/Simulate_Genotypes/4PopSplit/S1/C1/genos-test_common"
gwas_prefix = args[2] # Prefix to Test plink files
#gwas_prefix = "output/Simulate_Genotypes/4PopSplit/S1/C1/genos-gwas_common"
tvec_file = args[3] # Path to test vectors
#tvec_file = "output/Calculate_TGWAS/4PopSplit/S1/C1/Tvec.txt"
out_prefix = args[4] # Path to output directory
#out_prefix = "output/Calculate_TGWAS/4PopSplit/S1/C1/"

####################
## Functions #######
####################

# Compute G %*% t(X) %*% T
compute_b <- function(path_to_test, path_to_gwas, path_to_testvec, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  cmd_XT <- paste("sh code/Calculate_FGr/compute_XT.sh", path_to_test, path_to_testvec, outfile_XT, sep = " ")
  system(cmd_XT)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "xt_temp.Tvec.glm.linear"))
  count_plink <- fread(paste0(outpath, "xt_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  betas_plink_norm <- beta_plink$BETA * length_mc_genos

  # Compute GWAS genotype counts
  outfile_count <- paste0(outpath, "G_count")
  cmd_count <- paste("sh code/Calculate_FGr/compute_GWAScount.sh", path_to_gwas, outfile_count, sep = " ")
  system(cmd_count)

  # Calculate length of mean centered genotypes from counts
  count_plink <- fread(paste0(outpath, "G_count.gcount"))
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)
  length_mc_genos <- length_mc_genos * (1/(m-1))

  #  Re-write .linear file with correct betas
  beta_plink$BETA <- betas_plink_norm * (1/length_mc_genos)
  beta_reformat <- beta_plink %>% dplyr::select(ID, A1, BETA)
  fwrite(beta_reformat, paste0(outpath, "xt_temp.Tvec.glm.linear"), sep = "\t")

  # Compute b
  outfile_b <- paste0(outpath, "b")
  cmd_b <- paste("sh code/Calculate_FGr/GWAS_score.sh", path_to_gwas, paste0(outpath, "xt_temp.Tvec.glm.linear"), outfile_b, sep = " ")
  system(cmd_b)

  # Read in and return b
  b = fread(paste0(outpath, "b.sscore"))
  b = as.matrix(b$BETA_SUM)

  # Remove temporary files
  fn <- paste0(outpath, "xt_temp*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)
  fn <- paste0(outpath, "b*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)
  fn <- paste0(outpath, "G_count*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)

  return(b)
}

#####################
##     Main       ###
#####################

# Gather parameters
gwasID <- fread(paste0(gwas_prefix, ".psam"))
colnames(gwasID) <- c("FID", "IID",  "pheno")
m <- nrow(gwasID)
testID <- fread(paste0(test_prefix, ".psam"))
colnames(testID) <- c("FID", "IID",  "pheno")
n <- nrow(testID)

# Compute b
b = compute_b(path_to_test = test_prefix, path_to_gwas = gwas_prefix, path_to_testvec = tvec_file, outpath = out_prefix)
b = as.data.frame(scale(b))
colnames(b) <- "Tm"

fwrite(b, paste0(out_prefix, "Tm.txt"), row.names = F, col.names = T, quote = F, sep = "\t")




