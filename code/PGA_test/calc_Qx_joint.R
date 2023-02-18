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
tvec_file = args[3] # test vect
pops_file = args[4]
num = as.numeric(args[5]) # number of times to resapme
out_pre = args[6] # output prefix
true_file = args[7]
pc_list = args[8]
pc_list <- as.numeric(strsplit(pc_list,"-")[[1]])
print(pc_list)
cov_file = args[9]


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
  Va <- 4 * sum((effect_size)^2 * freq * (1 - freq))
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

# Function to calculate Q
calc_q <- function(sscore, Va, tvec) {

  tvec <- as.matrix(tvec)
  numerator <- t(tvec) %*% sscore
  qhat <- (1/sqrt(N * Va)) * (numerator)

  return(qhat)
}


# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(X, betas, Va, tvec) {

  # Flip effect sizes
  betas <- flip(betas)

  # Compute PGS
  prs <- pgs(X, betas)

  ## Calc Qx - Test
  q <- t(calc_q(prs, Va, tvec))

  return(q)
}


# Load Test vector
dfTvec <- fread(tvec_file)
Tvec <- dfTvec$Tvec - mean(dfTvec$Tvec)

# count number of TP individuals
N <- length(Tvec)

# Compute inverse longitude

## Load in covariance matrix
cov_mat <- fread(cov_file)
cov_mat <- apply(cov_mat, 2, as.numeric)

## Do eigen decomposition on covariance matrix
myE <-eigen(cov_mat)

## Compute inverse covariance matrix
n <- ncol(cov_mat) - 1
FXXinv <- myE$vectors[,1:n] %*% diag(1/sqrt(myE$values[1:n])) %*% t(myE$vectors[,1:n])

## Compute sd
gamma <- sd(FXXinv %*% Tvec)

## Rescale test vec
tvec_scaled <- Tvec / gamma

# Use rescaled test vec
tvec_scaled <- FXXinv %*% tvec_scaled


# Wrapper function to calculate Qx and empirical p values
main <- function(type, snps) {

  if(type == "TRUE") {
    # Load effect sizes
    beta_file <- true_file
    betas <- fread(beta_file, header=FALSE)
    colnames(betas) <- c("ID", "A1", "joint")
    betas$marginal <- betas$joint
    print(head(betas))
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

  # Mean center genotypes
  X <- scale(X, scale = FALSE)

  # Calc PGS
  sscore_marginal <- pgs(X, betas$marginal)
  sscore_joint <- pgs(X, betas$joint)

  # Calc Q
  q_marginal <- t(calc_q(sscore_marginal, Va_marginal, tvec_scaled))
  q_joint <- t(calc_q(sscore_joint, Va_joint, tvec_scaled))

  ## Calc Qx - Test
  #qx_marginal <- t(calc_Qx(sscore_marginal, tvec, Va_marginal, lambda_T))
  #qx_joint <- t(calc_Qx(sscore_joint, tvec, Va_joint, lambda_T))

  # Generate Empirical null marginal
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(X,  betas$marginal, Va_marginal, tvec_scaled)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_marginal <- length(all_strat[abs(all_strat) > abs(q_marginal[1,1])]) /length(all_strat)

  # Calculate p-value from chi-square marginal
  p_strat_marginal <- 2 * pnorm(abs(q_marginal[1,1]), lower.tail=FALSE)

  # Generate Empirical null joint
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(X, betas$joint, Va_joint, tvec_scaled)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_joint <- length(all_strat[abs(all_strat) > abs(q_joint[1,1])])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat_joint <- 2 * pnorm(abs(q_joint[1,1]), lower.tail=FALSE)

  # Concatenate output (Qx, p_strat)
  out <- c(q_marginal, p_strat_marginal , p_strat_en_marginal, q_joint, p_strat_joint, p_strat_en_joint)

  return(out)

}

# Run all types of PGS

### Asecertained
out_nc <- as.data.frame(matrix(NA, nrow = 4+length(pc_list), ncol =6))
out_nc[1, ] <- main(type = "TRUE")
out_nc[2, ] <- main(type = "", snps  = "nc")
out_nc[3, ] <- main(type = "-Tm", snps = "nc")
out_nc[4, ] <- main(type = "-ID", snps = "nc")
for (i in 1:length(pc_list)) {
  out_nc[(4+i), ]  <- main(type = paste0("-", pc_list[i]), snps = "nc")
}
pcs <- paste0("nc-PC", pc_list)
out_nc$type <- c("true","nc-uncorrected", "nc-Tm", "nc-ID", pcs)

### Causal
out_c <- as.data.frame(matrix(NA, nrow = 4+length(pc_list), ncol =6))
out_c[1, ] <- main(type = "TRUE")
out_c[2, ] <- main(type = "", snps  = "c")
out_c[3, ] <- main(type = "-Tm", snps = "c")
out_c[4, ] <- main(type = "-ID", snps = "c")
for (i in 1:length(pc_list)) {
  out_c[(4+i), ]  <- main(type = paste0("-",pc_list[i]), snps = "c")
}
pcs <- paste0("c-PC", pc_list)
out_c$type <- c("true","c-uncorrected", "c-Tm", "c-ID", pcs)

### Compute  bias
tq <- out_c[1,1]
if (is.na(tq)) {
   tq <- 0
}
out <- rbind(out_nc, out_c)
colnames(out) <- c("q_marginal", "P.Chi_marginal", "P.EN_marginal", "q_joint", "P.Chi_joint", "P.EN_joint", "type")
out$Bias_marginal <- out$q_marginal - tq
out$Bias_joint <- out$q_joint - tq


# Save output
print(out)
fwrite(out, out_pre,row.names=T,quote=F,sep="\t", col.names = T)


