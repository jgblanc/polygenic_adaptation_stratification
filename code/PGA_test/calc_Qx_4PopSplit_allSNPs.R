# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)!=9){stop("Rscript calc_TGWAS.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

gwas_prefix = args[1] # causal betas
geno_prefix = args[2] # Prefix to pilnk files
lambdaT_file = args[3]
tvec_file = args[4] # test vect
pops_file = args[5]
num = as.numeric(args[6]) # number of times to resapme
size = as.numeric(args[7]) # num individuals in the test panel
out_pre = args[8] # output prefix
true_file = args[9]


# Function to calculate Va
calc_Va <- function(afreq, es) {

  # Compute Va
  Va <- 2 * sum((es)^2 * afreq * (1 - afreq))

  return(Va)
}

# Function to calculate Qx
calc_Qx <- function(sscore, tvec, Va, lambda_T) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% sscore
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va*lambda_T)

  return(Qx_strat)
}

# Function to calculate Q
calc_q <- function(sscore, tvec) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% sscore

  return(Ztest)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(type, snps, tvec, Va, lambda_T) {

  if(type == "TRUE") {
    # Load effect sizes
    beta_file <- paste0(gwas_prefix, "true", ".betas")
    betas <- fread(beta_file)
    colnames(betas) <- c("ID", "A1", "BETA_Strat")
  } else {
    # Load effect sizes
    beta_file <- paste0(gwas_prefix, type,".", snps, ".betas")
    betas <- fread(beta_file)
    colnames(betas) <- c("ID", "A1", "BETA_Strat")
  }

  # Flip effect sizes
  betas$BETA_Strat <- flip(betas$BETA_Strat)

  # Save fliped effect sizes in temp file
  fwrite(betas, paste0(gwas_prefix,".tmp.betas"), row.names=F,quote=F,sep="\t", col.names = T)

  # Compute PGS using Plink
  prs_outfile <- paste0(gwas_prefix, ".tmp.prs")
  plink2_cmd <- paste("sh code/PGA_test/Test_score.sh", geno_prefix, paste0(gwas_prefix,".tmp.betas"), prs_outfile, sep = " ")
  system(plink2_cmd)
  sscore <- fread(paste0(prs_outfile, ".sscore")) %>% select(BETA_Strat_SUM)
  prs <- as.matrix(sscore)

  # Calculate Qx
  Qx <- t(calc_Qx(prs, tvec, Va, lambda_T))

  return(Qx)
}

# Load Lambda_T
lambda_T <- fread(lambdaT_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Load Test vector
std.tvec <- fread(tvec_file)
tvec <- std.tvec$Tvec / (nrow(std.tvec) - 1 )


# Wrapper function to calculate Qx and empirical p values
main <- function(type, snps) {

  if(type == "TRUE") {
    # Load effect sizes
    beta_file <- true_file
    betas <- fread(beta_file, header=FALSE)
    colnames(betas) <- c("ID", "A1", "BETA_Strat")
    beta_file <- paste0(gwas_prefix, "true", ".betas")
    fwrite(betas,beta_file,row.names=F,quote=F,sep="\t", col.names = T)
  } else {
    # Load effect sizes
    beta_file <- paste0(gwas_prefix, type,".", snps, ".betas")
    betas <- fread(beta_file)
    colnames(betas) <- c("ID", "A1", "BETA_Strat")
  }

  # Compute PGS using Plink
  prs_outfile <- paste0(gwas_prefix, ".prs")
  plink2_cmd <- paste("sh code/PGA_test/Test_score.sh", geno_prefix, beta_file, prs_outfile, sep = " ")
  system(plink2_cmd)
  sscore <- fread(paste0(prs_outfile, ".sscore")) %>% select(BETA_Strat_SUM)
  sscore <- as.matrix(sscore)

  # Compute Va
  freq <- fread(paste0(geno_prefix, ".afreq"))
  freq <- freq$ALT_FREQS
  Va <- calc_Va(afreq = freq, es = betas$BETA_Strat)

  ## Calc Q
  q <- t(calc_q(sscore, tvec))

  ## Calc Qx - Test
  qx <- t(calc_Qx(sscore, tvec, Va, lambda_T))

  # Generate Empirical null
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(type, snps, tvec, Va, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat <- pchisq(qx[1,1], df=1, lower.tail=FALSE)

  # Concatenate output (Qx, p_strat)
  out <- c(q, qx, p_strat , p_strat_en)

  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 7, ncol =4)
out[1, ] <- main(type = "", snps  = "nc")
out[2, ] <- main(type = "-Tm", snps = "nc")
out[3, ] <- main(type = "-ID", snps = "nc")
out[4, ] <- main(type = "", snps  = "c")
out[5, ] <- main(type = "-Tm", snps = "c")
out[6, ] <- main(type = "-ID", snps = "c")
out[7, ] <- main(type = "TRUE")
out <- as.data.frame(out)
out$Bias <- NA
tq <- out[7,1]
out[1,5] <- out[1,1] - tq
out[2,5] <- out[2,1] - tq
out[3,5] <- out[3,1] - tq
out[4,5] <- out[4,1] - tq
out[5,5] <- out[5,1] - tq
out[6,5] <- out[6,1] - tq
out[7,5] <- out[7,1] - tq

# Save output
colnames(out) <- c("Qx", "P.Chi", "P.EN", "Bias")
rownames(out) <- c("nc-uncorrected", "nc-Tm", "nc-ID", "c-uncorrected", "c-Tm", "c-ID", "true")
print(out)
fwrite(out, out_pre,row.names=T,quote=F,sep="\t", col.names = T)


