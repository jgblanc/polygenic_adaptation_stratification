## This script calculates Qx for all 3 types of GWAS (including and not including Tm)

## Requires: pgs (c, c.p, and n.c +/- Tm), va, va-Tm, lambda_T

args=commandArgs(TRUE)

if(length(args)<8){stop("Rscript calc_Tm.R <c.sscore> <c.p.sscore> <n.c.sscore> <lambda_T.txt> <Va.txt>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

c_file = args[1] # causal pgs
cp_file = args[2] # causal p-value pgs
nc_file = args[3] # clumped pgs
lambda_T_file = args[4] # lambda T
Va_file = args[5]
true_file = args[6] # True genetic value
tvec_file =args[7]
out_file = args[8] # outfile prefix

# Function to standardize PGS
stand_PGS <- function(prs_file, gv_file) {

  # Load PRS
  prs <- fread(prs_file)
  colnames(prs) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "RANDOM", "STRAT")

  # Load True GV
  gvalue <- fread(gv_file)
  colnames(gvalue) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "GV")

  # Join dataframes by IID
  df <- suppressMessages(full_join(prs,gvalue, by ="#IID")) %>% select("#IID", "RANDOM", "STRAT", "GV")

  # Standardize
  mprs.adj = df%>%
    mutate(random.adjusted = RANDOM-GV,
           strat.adjusted = STRAT-GV) %>%
    ungroup() %>% select("random.adjusted", "strat.adjusted")

  return(mprs.adj)
}

# Function to calculate Qx
calc_Qx <- function(mprs, tvec_file, Va, lambda_T) {
  #num_snps=200

  # Load Test vector
  std.tvec <- fread(tvec_file)
  n1 <- table(std.tvec)[1]
  n2 <- table(std.tvec)[2]
  std.tvec <- c(rep(1,(n1))/(n1), rep(-1,(n2))/(n2)) * (1/2)

  # Compute Qx Random
  Ztest <- t(std.tvec) %*% mprs$random.adjusted
  Qx_random <- (t(Ztest) %*% Ztest) / (Va[1]*lambda_T)

  # Compute Qx Strat
  Ztest <- t(std.tvec) %*% mprs$strat.adjusted
  Qx_strat <- (t(Ztest) %*% Ztest) / (Va[2]*lambda_T)

  # Calc p-values for chi-square df=1
  p_random <- pchisq(Qx_random, df=1, lower.tail=FALSE)
  p_strat <- pchisq(Qx_strat, df=1, lower.tail=FALSE)

  # Concat Results
  out <- c(Qx_random, Qx_strat, p_random, p_strat)
  return(out)
}

# Load Lambda_T
lambda_T <- fread(lambda_T_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Calculate Qx for standard GWAS
Va <- as.matrix(fread(Va_file), rownames=1)
Qx_mat <- matrix(NA, ncol = 4, nrow = 3)
colnames(Qx_mat) <- c("Qx_random", "Qx_strat", "p_random", "p_strat")
row.names(Qx_mat) <- rownames(Va)
Qx_mat[1,] <- calc_Qx(stand_PGS(c_file, true_file), tvec_file, Va[1,],lambda_T)
Qx_mat[2,] <- calc_Qx(stand_PGS(cp_file, true_file), tvec_file, Va[2,],lambda_T)
Qx_mat[3,] <- calc_Qx(stand_PGS(nc_file, true_file), tvec_file, Va[3,],lambda_T)

# Concat output
df <- as.data.frame(Qx_mat)
print(df)

# Save Va's
fwrite(df, out_file,row.names=T,quote=F,sep="\t", col.names = T)
