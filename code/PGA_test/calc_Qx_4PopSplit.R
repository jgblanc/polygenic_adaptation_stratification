# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)<16){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
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
pops_file = args[13]
num = as.numeric(args[14]) # number of times to resapme
out_pre = args[15] # output prefix
out_pgs = args[16]
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
  colnames(X) <- var.ids

  return(X)
}

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat_strat <- betas$BETA_Strat
  Z_strat <- X %*% Bhat_strat

  out <- cbind(Z_strat)
  colnames(out) <- c("STRAT")

  return(out)
}

# Function to standardize PGS
stand_PGS <- function(prs, gv_file) {

  # Load True GV
  gvalue <- fread(gv_file)
  colnames(gvalue) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "GV")

  # Join dataframes by IID
  df <- as.data.frame(cbind(prs, gvalue$GV))
  colnames(df) <- c("STRAT", "GV")

  # Standardize
  mprs.adj = df%>%
    mutate(strat.adjusted = STRAT-GV) %>%
    ungroup() %>% select("strat.adjusted")

  return(mprs.adj)
}

# Function to calculate Qx
calc_Qx <- function(mprs, tvec, Va, lambda_T) {

  # Compute Qx Strat
  Ztest <- t(tvec) %*% mprs$strat.adjusted
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va*lambda_T)

  return(Qx_strat)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(betas, tvec, Va, X, true_file, lambda_T) {

  # Flip effect sizes
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
Va_all <- as.matrix(fread(Va_file), rownames=1)
Va_Tm_all <- as.matrix(fread(Va_Tm_file), rownames=1)

# Load Test vector
std.tvec <- fread(tvec_file)
std.tvec <- std.tvec$V1
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
  pgs_stan <- stand_PGS(sscore, true_file)
  qx <- t(calc_Qx(pgs_stan, tvec, Va, lambda_T))

  # Calculate Ax
  mod <- lm(scale(tvec) ~ scale(pgs_stan$strat.adjusted))
  Ax <- coef(mod)[2]
  p_Ax <- summary(mod)$coefficients[2,4]

  # Generate Empirical null
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas, tvec, Va, X, true_file, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat <- pchisq(qx[1,1], df=1, lower.tail=FALSE)

  # Concatenate output (Qx_random, Qx_strat, p_random, p_strat)
  out <- c(qx, p_strat , p_strat_en, Ax, p_Ax)
  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 6, ncol =5)
out[1, ] <- main(c_file, Va_all[1,])
out[2, ] <- main(cp_file, Va_all[2,])
out[3, ] <- main(nc_file, Va_all[3,])
out[4, ] <- main(c_Tm_file, Va_Tm_all[1,])
out[5, ] <- main(cp_Tm_file, Va_Tm_all[2,])
out[6, ] <- main(nc_Tm_file, Va_Tm_all[3,])

# Save output
colnames(out) <- c("Qx", "P-Chi", "P-EN", "Ax", "P-Ax")
fwrite(out, out_pre,row.names=F,quote=F,sep="\t", col.names = T)

# Function to just output PGS
main2 <- function(beta_file) {

  # Load effect sizes
  betas <- fread(beta_file)

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)

  # Calc PGS
  sscore <- pgs(X, betas)

  # Calc Qx
  pgs_stan <- stand_PGS(sscore, true_file)

  return(pgs_stan$strat.adjusted)
}


# Output File with all the PGS
fam <- fread(paste0(geno_prefix, ".psam"))
fam <- fam[,1:2]
fam$c <- main2(c_file)
fam$c.p <- main2(cp_file)
fam$nc <- main2(nc_file)
fam$c_Tm <- main2(c_Tm_file)
fam$c.p_Tm <- main2(cp_Tm_file)
fam$nc_Tm <- main2(nc_Tm_file)
fam$Tvec <- tvec

fwrite(fam, out_pgs,row.names=F,quote=F,sep="\t", col.names = T)
