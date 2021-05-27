## This script calculates Qx for all 3 types of GWAS (including and not including Tm)

## Requires: pgs (c, c.p, and n.c +/- Tm), va, va-Tm, lambda_T
args=commandArgs(TRUE)

if(length(args)<12){stop("Rscript calc_Tm.R <c.sscore> <c.p.sscore> <n.c.sscore> <lambda_T.txt> <Va.txt>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

c_file = args[1] # causal pgs
cp_file = args[2] # causal p-value pgs
nc_file = args[3] # clumped pgs
c_Tm_file = args[4] # causal pgs (Tm)
cp_Tm_file = args[5] # causal p-value pgs (Tm)
nc_Tm_file = args[6] # clumped pgs (Tm)
lambda_T_file = args[7] # lambda T
Va_file = args[8]
Va_Tm_file = args[9]
true_file = args[10] # True genetic value
tvec_file =args[11]
resamples_file = args[12]
out_file = args[13] # outfile prefix

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
calc_Qx <- function(mprs, tvec_file, Va, lambda_T, resamples) {

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
  all_random <- dplyr::pull(resamples, 1)
  p_random <- length(all_random[all_random > Qx_random[1,1]])/length(all_random)
  all_strat <- dplyr::pull(resamples, 2)
  p_strat <- length(all_strat[all_strat > Qx_strat[1,1]])/length(all_strat)

  # Concat Results
  out <- c(Qx_random, Qx_strat, p_random, p_strat)
  return(out)
}

# Load Lambda_T
lambda_T <- fread(lambda_T_file)
lambda_T <- as.numeric(as.character(lambda_T[1,1]))

# Load flipped sign resamples
flip <- fread(resamples_file)

# Calculate Qx for standard GWAS
Va <- as.matrix(fread(Va_file), rownames=1)
Qx_mat <- matrix(NA, ncol = 4, nrow = 3)
colnames(Qx_mat) <- c("Qx_random", "Qx_strat", "p_random", "p_strat")
row.names(Qx_mat) <- rownames(Va)
Qx_mat[1,] <- calc_Qx(stand_PGS(c_file, true_file), tvec_file, Va[1,],lambda_T, flip[,1:2])
Qx_mat[2,] <- calc_Qx(stand_PGS(cp_file, true_file), tvec_file, Va[2,],lambda_T,flip[,3:4] )
Qx_mat[3,] <- calc_Qx(stand_PGS(nc_file, true_file), tvec_file, Va[3,],lambda_T, flip[,5:6])

# Calculate Qx for Tm GWAS
Va <- as.matrix(fread(Va_Tm_file), rownames=1)
Qx_mat_Tm <- matrix(NA, ncol = 4, nrow = 3)
colnames(Qx_mat_Tm) <- c("Qx_random", "Qx_strat", "p_random", "p_strat")
row.names(Qx_mat_Tm) <- c("Tm-c", "Tm-c.p", "Tm-nc")
Qx_mat_Tm[1,] <- calc_Qx(stand_PGS(c_Tm_file, true_file), tvec_file, Va[1,],lambda_T, flip[,7:8])
Qx_mat_Tm[2,] <- calc_Qx(stand_PGS(cp_Tm_file, true_file), tvec_file, Va[2,],lambda_T, flip[,9:10])
Qx_mat_Tm[3,] <- calc_Qx(stand_PGS(nc_Tm_file, true_file), tvec_file, Va[3,],lambda_T, flip[,11:12])


# Concat output
df <- as.data.frame(rbind(Qx_mat, Qx_mat_Tm))
fwrite(df, out_file,row.names=T,quote=F,sep="\t", col.names = T)



