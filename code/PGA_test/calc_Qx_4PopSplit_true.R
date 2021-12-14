# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

es_file = args[1] # true effect sizes
geno_prefix = args[2] # test panel genotypes
tvec_file = args[3] # test vector
lambdaT_file = args[4]
pops_file = args[5]
num = as.numeric(args[6]) # number of times to resapme
out_pre = args[7]


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

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat <- betas$BETA
  Z <- X %*% Bhat

  out <- cbind(Z)
  colnames(out) <- c("GV")

  return(out)
}

# Function to calculate Va
calc_Va <- function(geno_mat, es) {

  # Get allele frequency in test panel
  freq <- colMeans(geno_mat) /2

  # Pull out effect sizes only
  effect_size <- es$BETA

  # Compute Va
  Va <- 2 * sum((effect_size)^2 * freq * (1 - freq))
  return(Va)
}


# Function to calculate Qx
calc_Qx <- function(mprs, tvec, Va, lambda_T) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% mprs
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va*lambda_T)

  return(Qx_strat)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, lambda_T) {

  # Flip effect sizes
  betas$BETA <- flip(betas$BETA)

  # Calculate PGS
  prs <- pgs(X, betas)

  # Calculate Qx
  Qx <- t(calc_Qx(prs, tvec, Va, lambda_T))

  return(Qx)
}


# Load Test vector
std.tvec <- fread(tvec_file)
std.tvec <- std.tvec$V1
n1 <- table(std.tvec)[1]
n2 <- table(std.tvec)[2]
tvec <- c(rep(1,(n1))/(n1), rep(-1,(n2))/(n2)) * (1/2)

# Load Lambda_T
lambda_T <- fread(lambdaT_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))


# Wrapper function to calculate Qx and empirical p values
main <- function(beta_file) {

  # Load effect sizes
  betas <- fread(beta_file)
  colnames(betas) <- c("ID", "A1", "BETA")

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)


  # Calc PGS
  sscore <- pgs(X, betas)

  # Calc Va
  Va <- calc_Va(X, betas)

  # Calc Qx
  qx <- t(calc_Qx(sscore, tvec, Va, lambda_T))

  # Generate Empirical null
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas, tvec, Va, X, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat <- pchisq(qx[1,1], df=1, lower.tail=FALSE)

  # Concatenate output (Qx_random, Qx_strat, p_random, p_strat)
  out <- c(qx, p_strat , p_strat_en)
  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 1, ncol =3)
out[1, ] <- main(es_file)
print(out)

# Save output
colnames(out) <- c("Qx", "P-Chi", "P-EN")
fwrite(out, out_pre,row.names=F,quote=F,sep="\t", col.names = T)


