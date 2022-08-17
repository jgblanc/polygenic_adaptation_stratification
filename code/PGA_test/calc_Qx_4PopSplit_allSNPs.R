# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)!=8){stop("Rscript calc_TGWAS.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
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
  effect_size <- es$BETA

  # Compute Va
  Va <- 2 * sum((effect_size)^2 * freq * (1 - freq))
  return(Va)
}


# Function to compute PGS
pgs <- function(X, betas) {

  Bhat_strat <- betas$BETA_Strat
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

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, lambda_T) {

  # Flip effect sizes
  betas$BETA_Strat <- flip(betas$BETA_Strat)

  # Calculate PGS
  prs <- pgs(X, betas)

  # Calculate Qx
  Qx <- t(calc_Qx(prs, tvec, Va, lambda_T))

  return(Qx)
}

# Load Lambda_T
lambda_T <- fread(lambdaT_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Load Test vector
std.tvec <- fread(tvec_file)
tvec <- std.tvec$Tvec


# Wrapper function to calculate Qx and empirical p values
main <- function(type) {

  Va <- rep(0,200)
  ssMat <- matrix(NA, nrow = 200,ncol=size)
  for (i in 1:200) {

    # Load effect sizes from CHR i-1
    beta_file <- paste0(gwas_prefix, "_", (i-1), type, ".betas")
    print(beta_file)
    betas <- fread(beta_file)
    head(betas)
    colnames(betas) <- c("ID", "A1", "BETA_Strat")

    # Load Genotypes
    X <- read_genos(geno_prefix, betas)

    # Calculate Va
    Va[i] <- calc_Va(X, betas)

    # Calc PGS
    ssMat[i,] <- pgs(X, betas)

  }

  Va <- sum(Va)
  print(Va)
  sscore <- colSums(ssMat)
  print(sscore)

  ## Calc Qx - Test
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

  # Concatenate output (Qx, p_strat)
  out <- c(qx, p_strat , p_strat_en)
  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 4, ncol =3)
out[1, ] <- main(type = "")
out[2, ] <- main(type = "-Tm")
out[3, ] <- main(type = "-ID")
print(out)

# Save output
colnames(out) <- c("Qx", "P.Chi", "P.EN")
rownames(out) <- c("uncorrected", "Tm", "ID")
print(out)
fwrite(out, out_pre,row.names=F,quote=F,sep="\t", col.names = T)


