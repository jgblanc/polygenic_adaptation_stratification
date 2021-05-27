# This script takes the Va, .sscore, and Fst for each resamples replicate and calculates Qx


args=commandArgs(TRUE)

if(length(args)<16){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

c_file = args[1] # causal pgs
cp_file = args[2] # causal p-value pgs
nc_file = args[3] # clumped pgs
c_Tm_file = args[4] # causal pgs Tm
cp_Tm_file = args[5] # causal p-value pgs Tm
nc_Tm_file = args[6] # clumped pgs Tm
lambda_file = args[7] # Lambda T
va_c_file = args[8]
va_cp_file = args[9]
va_nc_file = args[10]
va_c_Tm_file = args[11]
va_cp_Tm_file = args[12]
va_nc_Tm_file = args[13]
true_file = args[14]
tvec_file = args[15]
out_pre = args[16]

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
calc_Qx <- function(mprs, tvec_file, Va_file, lambda_T) {

  # Load Va
  Va <- fread(Va_file)
  Va <- as.numeric(Va[1,])

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

  # Concat Results
  out <- c(Qx_random, Qx_strat)
  return(out)
}

# Load Lambda_T
lambda_T <- fread(lambda_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Calculate Qx
qx_c <- t(calc_Qx(stand_PGS(c_file, true_file), tvec_file, va_c_file, lambda_T))
qx_cp <- t(calc_Qx(stand_PGS(cp_file, true_file), tvec_file, va_cp_file, lambda_T))
qx_nc <- t(calc_Qx(stand_PGS(nc_file, true_file), tvec_file, va_nc_file, lambda_T))
qx_c_Tm <- t(calc_Qx(stand_PGS(c_Tm_file, true_file), tvec_file, va_c_Tm_file, lambda_T))
qx_cp_Tm <- t(calc_Qx(stand_PGS(cp_Tm_file, true_file), tvec_file, va_cp_Tm_file, lambda_T))
qx_nc_Tm <- t(calc_Qx(stand_PGS(nc_Tm_file, true_file), tvec_file, va_nc_Tm_file, lambda_T))

# Output all qx
df <- cbind(qx_c, qx_cp, qx_nc, qx_c_Tm, qx_cp_Tm, qx_nc_Tm)
fwrite(df, out_pre,row.names=F,quote=F,sep="\t", col.names = F)

