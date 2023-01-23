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


# Function to read in genotype matrix for a set of variants
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Function to calculate Va
calc_Va <- function(geno_mat, es) {

  # Get allele frequency in test panel
  freq <- colMeans(geno_mat) /2

  # Pull out effect sizes only
  effect_size <- es

  # Compute Va
  Va <- 2 * sum((effect_size)^2 * freq * (1 - freq))
  return(Va)
}

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat_strat <- betas
  Z_strat <- X %*% Bhat_strat

  out <- cbind(Z_strat)
  colnames(out) <- c("STRAT")

  return(out)
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
en <- function(betas, tvec, Va, lambda_T) {

  # Flip effect sizes
  betas <- flip(betas)

  # Compute PGS
  prs <- pgs(X, betas)


  ## Calc Qx - Test
  qx <- t(calc_Qx(prs, tvec, Va, lambda_T))

  return(qx)
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
    colnames(betas) <- c("ID", "A1", "joint")
    betas$marginal <- betas$joint
  } else {
    # Load effect sizes
    beta_file <- paste0(gwas_prefix, type,".", snps, ".betas.joint")
    betas <- fread(beta_file, header = TRUE)
  }

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)

  # Calculate Va
  Va_marginal <- calc_Va(X, betas$marginal)
  Va_joint <- calc_Va(X, betas$joint)

  # Calc PGS
  sscore_marginal <- pgs(X, betas$marginal)
  sscore_joint <- pgs(X, betas$joint)

  ## Calc Q
  q_marginal <- t(calc_q(sscore_marginal, tvec))
  q_joint <- t(calc_q(sscore_joint, tvec))

  ## Calc Qx - Test
  qx_marginal <- t(calc_Qx(sscore_marginal, tvec, Va, lambda_T))
  qx_joint <- t(calc_Qx(sscore_joint, tvec, Va, lambda_T))

  # Generate Empirical null marginal
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas$marginal, tvec, Va_marginal, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_marginal <- length(all_strat[all_strat > qx_marginal[1,1]])/length(all_strat)

  # Calculate p-value from chi-square marginal
  p_strat_marginal <- pchisq(qx_marginal[1,1], df=1, lower.tail=FALSE)

  # Generate Empirical null joint
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas$joint, tvec, Va_joint, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_joint <- length(all_strat[all_strat > qx_joint[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat_joint <- pchisq(qx_joint[1,1], df=1, lower.tail=FALSE)

  # Concatenate output (Qx, p_strat)
  out <- c(q_marginal, qx_marginal, p_strat_marginal , p_strat_en_marginal, q_joint, qx_joint, p_strat_joint, p_strat_en_joint)

  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 7, ncol =8)
out[1, ] <- main(type = "", snps  = "nc")
out[2, ] <- main(type = "-Tm", snps = "nc")
out[3, ] <- main(type = "-ID", snps = "nc")
out[4, ] <- main(type = "", snps  = "c")
out[5, ] <- main(type = "-Tm", snps = "c")
out[6, ] <- main(type = "-ID", snps = "c")
out[7, ] <- main(type = "TRUE")
out <- as.data.frame(out)
out$Bias_marginal <- NA
out$Bias_joint <- NA
tq <- out[7,1]

## Marginal bias
out[1,9] <- out[1,1] - tq
out[2,9] <- out[2,1] - tq
out[3,9] <- out[3,1] - tq
out[4,9] <- out[4,1] - tq
out[5,9] <- out[5,1] - tq
out[6,9] <- out[6,1] - tq
out[7,9] <- out[7,1] - tq

## Joint bias
out[1,10] <- out[1,5] - tq
out[2,10] <- out[2,5] - tq
out[3,10] <- out[3,5] - tq
out[4,10] <- out[4,5] - tq
out[5,10] <- out[5,5] - tq
out[6,10] <- out[6,5] - tq
out[7,10] <- out[7,5] - tq

# Save output
colnames(out) <- c("q_marginal","qx_marginal", "P.Chi_marginal", "P.EN_marginal", "q_joint", "qx_joint", "P.Chi_joint", "P.EN_joint", "Bias_margina", "Bias_joint")
rownames(out) <- c("nc-uncorrected", "nc-Tm", "nc-ID", "c-uncorrected", "c-Tm", "c-ID", "true")
print(out)
fwrite(out, out_pre,row.names=T,quote=F,sep="\t", col.names = T)


