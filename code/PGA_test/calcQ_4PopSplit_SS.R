# This script takes a setof effect sizes and test vector and computes Q

args=commandArgs(TRUE)

if(length(args)!=8){stop("Rscript calcQ_4PopSplit.R <prefix to effect sizes> <test panel genotypes prefix> <test vector>
                         <number of times to resample in empirical null> <output prefix> <true effect sizes> <list of PCs>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

geno_test_prefix = args[1] # Prefix to pilnk files
tvec_file = args[2] # test vect
pops_file = args[3]
num = as.numeric(args[4]) # number of times to resapme
out_pre = args[5] # output prefix
true_file = args[6]
out_pgs <- args[7]
betas_file = args[8]
outpath = args[9]
geno_prefix = geno_test_prefix


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

## Compute r

# Load Test vector
dfTvec <- fread(tvec_file)
Tvec <- dfTvec$Tvec - mean(dfTvec$Tvec)

# count number of TP individuals
N <- length(Tvec)

# Rescale test vec to have variance 1
tvec_scaled <- scale(Tvec)

# Compute t(X)T using plink
outfile_XT <- paste0(outpath, "xt")
cmd_XT <- paste("sh code/Calculate_FGr/compute_XT.sh", geno_test_prefix, tvec_file, outfile_XT, sep = " ")
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
r <- beta_reformat$BETA / (N -1)

## Compute corrected effect sizes
betas <- fread(betas_file)
colnames(betas)[1] <- "CHROM"

# Function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  gwas_df = gwas_df[A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P, SE)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

# Flip affect to "T"
betas <- flip_effect(betas, "BETA")

# Regress r out of betas
mod <- lm(betas$BETA ~ r)
betas$corrected <- mod$residuals

# Convert to Chi-sq
#r_se <- sd(r) / sqrt(length(r))
#Beta_r <- mod$coefficients[2]
#seBr <-   (betas$SE)^2 + (Beta_r^2 * (beta_plink$SE^2))
#chsq2 <- (mod$residuals)^2 /(seBr^2)
#pvals <- 1- pchisq(chsq2, 1)


# Clump snps - random
#betas <- betas %>% mutate(ID2 = ID) %>%
#  separate(ID2, c("chr", "pos", "a1", "a2"), "_")%>% group_by(chr) %>%
#  slice_min(P, with_ties = F) %>% ungroup() %>%
#  select(CHROM, POS, ID, A1, BETA, P, corrected)

#betas <- betas %>% mutate(ID2 = ID) %>%
#  separate(ID2, c("chr", "pos", "a1", "a2"), "_")%>% group_by(chr) %>%
#  sample_n(1) %>% ungroup() %>%
#  select(CHROM, POS, ID, A1, BETA, P, corrected)

# Get causal sites
ts <- fread(true_file)
colnames(ts) <- c("ID", "A1", "TS")
betas_ts <- inner_join(ts, betas)
betas_ts$BETA <- betas_ts$TS
betas_ts$corrected <- betas_ts$TS

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
  numerator <-  (1/(N-1)) * t(tvec) %*% sscore
  #qhat <- (sqrt(N / Va)) * (numerator)
  qhat <- numerator

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

# Rescale test vec to have variance 1
tvec_scaled <- scale(Tvec)

# Wrapper function to calculate Qx and empirical p values
main <- function(betas) {

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)

  # Calculate Va
  Va_uc <- calc_Va(X, betas$BETA)
  Va_c <- calc_Va(X, betas$corrected)

  # Mean center genotypes
  X <- scale(X, scale = FALSE)

  # Calc PGS
  sscore_uc <- pgs(X, betas$BETA)
  sscore_c <- pgs(X, betas$corrected)

  # Calc Q
  q_uc <- t(calc_q(sscore_uc, Va_uc, tvec_scaled))
  q_c <- t(calc_q(sscore_c, Va_c, tvec_scaled))

  # Generate Empirical null marginal
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(X,  betas$BETA, Va_uc, tvec_scaled)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_uc <- length(all_strat[abs(all_strat) > abs(q_uc[1,1])]) /length(all_strat)

  # Generate Empirical null joint
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(X, betas$corrected, Va_c, tvec_scaled)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_c <- length(all_strat[abs(all_strat) > abs(q_c[1,1])])/length(all_strat)

  # Concatenate output (Qx, p_strat)
  out <- c(q_uc, p_strat_en_uc, q_c, p_strat_en_c)

  return(out)

}

# Run all types of PGS

### Causal
out <- as.data.frame(matrix(NA, nrow = 2, ncol =4))
out[1, ] <- main(betas_ts)
out[2, ] <- main(betas)
out$type <- c("true","nc-uncorrected")


### Compute  bias
tq <- out[1,1]
if (is.na(tq)) {
   tq <- 0
}
out <- rbind(out_nc, out_c)
colnames(out) <- c("q_uc", "P.EN_uc", "q_c", "P.EN_c", "type")
out$Bias_marginal <- out$q_marginal - tq
out$Bias_joint <- out$q_joint - tq
out$Bias_uc <- out$q_uc - tq
out$Bias_c <- out$q_c - tq


# Save output
print(out)
fwrite(out, out_pre,row.names=T,quote=F,sep="\t", col.names = T)


