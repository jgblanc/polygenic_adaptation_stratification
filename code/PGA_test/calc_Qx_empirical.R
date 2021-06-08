# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)<14){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

c_file = args[1] # causal betas
cp_file = args[2] # causal p-value betas
nc_file = args[3] # clumped betas
c_Tm_file = args[4] # causal betas
cp_Tm_file = args[5] # causal p-value betas
nc_Tm_file = args[6] # clumped betas
geno_prefix = args[7] # Prefix to pilnk files
lambdaT_file = args[8]
Va_file = args[9] # Va
Va_Tm_file = args[10] # Va Tm
true_file = args[11] # true PGS
tvec_file = args[12] # test vect
num = as.numeric(args[13]) # number of times to resapme
out_pre = args[14] # output prefix
#out_pre = "~/polygenic_adaptation_stratification/output/Empirical_Null/4PopSplit/E1/C1/h2-0/env-0.0/geno-gwas_"


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

  return(X)
}

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat_random <- betas$BETA_Random
  Z_random <- X %*% Bhat_random

  Bhat_strat <- betas$BETA_Strat
  Z_strat <- X %*% Bhat_strat

  out <- cbind(Z_random, Z_strat)
  colnames(out) <- c("RANDOM", "STRAT")

  return(out)
}

# Function to standardize PGS
stand_PGS <- function(prs, gv_file) {

  # Load True GV
  gvalue <- fread(gv_file)
  colnames(gvalue) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "GV")

  # Join dataframes by IID
  df <- as.data.frame(cbind(prs, gvalue$GV))
  colnames(df) <- c("RANDOM", "STRAT", "GV")

  # Standardize
  mprs.adj = df%>%
    mutate(random.adjusted = RANDOM-GV,
           strat.adjusted = STRAT-GV) %>%
    ungroup() %>% select("random.adjusted", "strat.adjusted")

  return(mprs.adj)
}

# Function to calculate Qx
calc_Qx <- function(mprs, tvec, Va, lambda_T) {
  std.tvec <- tvec

  # Compute Qx Random
  Ztest <- t(std.tvec) %*% mprs$random.adjusted
  Qx_random <- (t(Ztest) %*% Ztest) / (Va[1]*lambda_T)

  # Compute Qx Strat
  Ztest <- t(std.tvec) %*% mprs$strat.adjusted
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va[2]*lambda_T)

  # Concat Results
  out <- c(Qx_random, Qx_strat)
  return(out)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, true_file, lambda_T) {

  # Flip effect sizes
  betas$BETA_Random <- flip(betas$BETA_Random)
  betas$BETA_Strat <- flip(betas$BETA_Strat)

  # Calculate PGS
  prs <- pgs(X, betas)

  # Calculate Qx
  Qx <- t(calc_Qx(stand_PGS(prs, true_file), tvec, Va, lambda_T))

  return(Qx)
}

# Load Lambda_T
lambda_T <- fread(lambdaT_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Load Va
Va <- as.matrix(fread(Va_file), rownames=1)
Va_Tm <- as.matrix(fread(Va_Tm_file), rownames=1)

# Load Test vector
std.tvec <- fread(tvec_file)
n1 <- table(std.tvec)[1]
n2 <- table(std.tvec)[2]
tvec <- c(rep(1,(n1))/(n1), rep(-1,(n2))/(n2)) * (1/2)


# Wrapper function to calculate Qx and empirical p values
main <- function(beta_file, Va) {

  # Load effect sizes
  betas <- fread(beta_file)
  #colnames(betas) <- c("ID", "A1", "BETA_Random", "BETA_Strat")  

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)
 
  # Calc PGS
  sscore <- pgs(X, betas)

  # Calc Qx
  qx <- t(calc_Qx(stand_PGS(sscore, true_file), tvec, Va, lambda_T))

  # Generate Empirical null
  redraws <- matrix(0, ncol = 2, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas, tvec, Va, X, true_file, lambda_T)
  }

  # Calculate p-values
  all_random <- redraws[,1]
  p_random <- length(all_random[all_random > qx[1,1]])/length(all_random)
  all_strat <- redraws[,2]
  p_strat <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Concatenate output (Qx_random, Qx_strat, p_random, p_strat)
  out <- c(qx, p_random, p_strat)
  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 6, ncol =4)
out[1, ] <- main(c_file, Va[1,])
out[2, ] <- main(cp_file, Va[2,])
out[3, ] <- main(nc_file, Va[3,])
out[4, ] <- main(c_Tm_file, Va_Tm[1,])
out[5, ] <- main(cp_Tm_file, Va_Tm[2,])
out[6, ] <- main(nc_Tm_file, Va_Tm[3,])

# Save output
fwrite(out, out_pre,row.names=F,quote=F,sep="\t", col.names = F)


