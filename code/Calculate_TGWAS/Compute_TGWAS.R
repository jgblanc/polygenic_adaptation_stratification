# This script uses conjugate gradient descent to solve for TGWAS

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript Compute_TGWAS.R <test_prefix> <gwas_prefix> <tvec_file> <out_prefix>")}

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
print(out_prefix)

####################
## Functions #######
####################

# Compute G %*% t(X) %*% T
compute_b <- function(path_to_test, path_to_gwas, path_to_testvec, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  cmd_XT <- paste("sh code/Calculate_TGWAS/compute_XT.sh", path_to_test, path_to_testvec, outfile_XT, sep = " ")
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

  #  Re-write .linear file with correct betas
  beta_plink$BETA <- betas_plink_norm
  beta_reformat <- beta_plink %>% dplyr::select(ID, A1, BETA)
  fwrite(beta_reformat, paste0(outpath, "xt_temp.Tvec.glm.linear"), sep = "\t")

  # Compute b
  outfile_b <- paste0(outpath, "b")
  cmd_b <- paste("sh code/Calculate_TGWAS/GWAS_score.sh", path_to_gwas, paste0(outpath, "xt_temp.Tvec.glm.linear"), outfile_b, sep = " ")
  system(cmd_b)

  # Read in and return b
  b = fread(paste0(outpath, "b.sscore"))
  b = as.matrix(b$BETA_SUM)
  return(b)
}

# Compute residual
compute_rk <- function(x, path_to_gwas, outpath, gwasID) {

  # Create pheno file for Xk
  gwasID$pheno <- x
  fwrite(gwasID, paste0(outpath,"xk.txt"), sep = "\t")

  # Use plink to comput t(G)x
  outfile_Gx <- paste0(outpath, "Gx_temp")
  path_to_pheno <- paste0(outpath,"xk.txt")
  cmd_Gx <- paste("sh code/Calculate_TGWAS/GWAS_assoc.sh", path_to_gwas, path_to_pheno, outfile_Gx, sep = " ")
  system(cmd_Gx)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "Gx_temp.pheno.glm.linear"))
  count_plink <- fread(paste0(outpath, "Gx_temp.gcount"))

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
  fwrite(beta_reformat, paste0(outpath, "Gx_temp.pheno.glm.linear"), sep = "\t")

  # Compute GGx
  outfile_ggx <- paste0(outpath, "GGx")
  cmd_ggx <- paste("sh code/Calculate_TGWAS/GWAS_score.sh", path_to_gwas, paste0(outpath, "Gx_temp.pheno.glm.linear"), outfile_ggx, sep = " ")
  system(cmd_ggx)

  # Compute rk and return
  GGx <- fread(paste0(outpath, "GGx.sscore"))
  b <- fread(paste0(outpath, "b.sscore"))
  rk <- GGx$BETA_SUM - b$BETA_SUM
  return(rk)
}

# Compute next step size (alpha)
compute_alpha <- function(r, p, path_to_gwas, outpath, gwasID) {

  # Create pheno file for pk
  gwasID$pheno <- p
  fwrite(gwasID, paste0(outpath,"pk.txt"), sep = "\t")

  # Compute t(G) %*% pk using plink
  outfile_Gp <- paste0(outpath, "Gp_temp")
  path_to_pheno <- paste0(outpath,"pk.txt")
  cmd_Gp <- paste("sh code/Calculate_TGWAS/GWAS_assoc.sh", path_to_gwas, path_to_pheno, outfile_Gp, sep = " ")
  system(cmd_Gp)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "Gp_temp.pheno.glm.linear"))
  count_plink <- fread(paste0(outpath, "Gp_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  betas_plink_norm <- beta_plink$BETA * length_mc_genos

  # Compute and return alpha
  alpha <- (t(r) %*% r) / sum(betas_plink_norm^2)
  return(alpha)
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

# Initialize X0
x0 = as.matrix(rep(0, m))
xk = x0
# Compute r0
rk = -b
# Assign initial direction
pk = -rk
# Initialize stop criteria and iteration counter
sc = 1
counter = 1
while (sc > 1e-5) {

  if (counter > 100) {
    break
  }

  print(paste0("Iteration: ", counter))

  # Compute rkrk
  rkrk = t(rk) %*% rk

  # Compute alpha
  alpha_k = compute_alpha(r = rk, p = pk, path_to_gwas = gwas_prefix, outpath = out_prefix, gwasID = gwasID)

  # Compute X_k+1
  xk_prev = xk
  xk = xk + (alpha_k[1,1] * pk)

  # Compute new residual
  rk = compute_rk(x = xk, path_to_gwas = gwas_prefix , outpath = out_prefix, gwasID = gwasID)

  # Compute Beta_k
  beta_k = t(rk) %*% rk / rkrk

  # Compute p_k+1
  pk = -rk + beta_k[1,1] * pk

  # Update counters
  counter = counter + 1
  sc = norm(xk - xk_prev, type = "2") / norm(xk_prev, type = "2")
  print(sc)

}

fwrite(xk, paste0(out_prefix, "TGWAS.txt"), row.names = F, col.names = F, quote = F, sep = "\t")




