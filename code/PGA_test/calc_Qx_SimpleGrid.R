# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null

args=commandArgs(TRUE)

if(length(args)!=17){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
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
c_ID_file = args[7] # causal betas
cp_ID_file = args[8] # causal p-value betas
nc_ID_file = args[9] # clumped betas
geno_prefix = args[10] # Prefix to pilnk files
lambdaT_file = args[11]
true_file = args[12] # true PGS
tvec_file = args[13] # test vect
pops_file = args[14]
num = as.numeric(args[15]) # number of times to resapme
out_pre = args[16] # output prefix
out_pgs = args[17]

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

# Function to clean PGS
stand_PGS <- function(prs, gv_file) {

  # Load True GV
  gvalue <- fread(gv_file)
  colnames(gvalue) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "GV")

  # Join dataframes by IID
  df <- as.data.frame(cbind(prs, gvalue$GV))
  colnames(df) <- c("STRAT", "GV")

  mprs.adj = df%>%
    mutate(strat.adjusted = STRAT) %>%
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

# Load Test vector
std.tvec <- fread(tvec_file)
tvec <- std.tvec$V1

# Make longitude test vector
pops <- fread(pops_file)
fam <-  fread(paste0(geno_prefix, ".psam"))
colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))
Tvec <- pop$Long
tvec_long <- (Tvec-mean(Tvec))

# Make latitude test vector
pops <- fread(pops_file)
fam <-  fread(paste0(geno_prefix, ".psam"))
colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))
Tvec <- pop$Lat
tvec_lat <- (Tvec-mean(Tvec))

# Wrapper function to calculate Qx and empirical p values
main <- function(beta_file) {

  # Load effect sizes
  betas <- fread(beta_file)

  # Load Genotypes
  X <- read_genos(geno_prefix, betas)

  # Calculate Va
  Va <- calc_Va(X, betas)

  # Calc PGS
  sscore <- pgs(X, betas)
  pgs_stan <- stand_PGS(sscore, true_file)

  ## Calc Qx - Test
  qx <- t(calc_Qx(pgs_stan, tvec, Va, lambda_T))

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


  ## Calc Qx - Lat
  qx <- t(calc_Qx(pgs_stan, tvec_lat, Va, lambda_T))

  # Generate Empirical null
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas, tvec, Va, X, true_file, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_lat <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat_lat <- pchisq(qx[1,1], df=1, lower.tail=FALSE)

  ## Calc Qx - Long
  qx <- t(calc_Qx(pgs_stan, tvec_long, Va, lambda_T))

  # Generate Empirical null
  redraws <- matrix(0, ncol = 1, nrow = num)
  for (i in 1:num){
    redraws[i,] <- en(betas, tvec, Va, X, true_file, lambda_T)
  }

  # Calculate empirical p-values
  all_strat <- redraws[,1]
  p_strat_en_long <- length(all_strat[all_strat > qx[1,1]])/length(all_strat)

  # Calculate p-value from chi-square
  p_strat_long <- pchisq(qx[1,1], df=1, lower.tail=FALSE)

  # Concatenate output (Qx_random, Qx_strat, p_random, p_strat)
  out <- c(qx, p_strat , p_strat_en, p_strat_lat, p_strat_en_lat, p_strat_long, p_strat_en_long)
  return(out)

}

# Run all types of PGS
out <- matrix(NA, nrow = 9, ncol =7)
out[1, ] <- main(c_file)
out[2, ] <- main(cp_file)
out[3, ] <- main(nc_file)
out[4, ] <- main(c_Tm_file)
out[5, ] <- main(cp_Tm_file)
out[6, ] <- main(nc_Tm_file)
out[7, ] <- main(c_ID_file)
out[8, ] <- main(cp_ID_file)
out[9, ] <- main(nc_ID_file)
print(out)

# Save output
colnames(out) <- c("Qx", "P.Chi", "P.EN", "P.Chi_Lat", "P.EN_Lat", "P.Chi_Long", "P.EN_Long")
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
fam$c_ID <- main2(c_ID_file)
fam$c.p_ID <- main2(cp_ID_file)
fam$nc_ID <- main2(nc_ID_file)
fam$Tvec <- tvec

fwrite(fam, out_pgs,row.names=F,quote=F,sep="\t", col.names = T)
