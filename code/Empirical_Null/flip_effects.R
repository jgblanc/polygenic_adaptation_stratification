# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null
# It also re-calculates Va for these new effect sizes

args=commandArgs(TRUE)

if(length(args)<9){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

c_file = args[1] # causal betas
cp_file = args[2] # causal p-value betas
nc_file = args[3] # clumped betas
c_Tm_file = args[4] # causal betas
cp_Tm_file = args[5] # causal p-value betas
nc_Tm_file = args[6] # clumped betas
freq_file = args[7]
num = args[8] # number of times to resample
out_pre = args[9] # output prefix
#out_pre = "~/polygenic_adaptation_stratification/output/Empirical_Null/4PopSplit/E1/C1/h2-0/env-0.0/geno-gwas_"

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to write new betas to a file
write_file <- function(betas_df, betas_random, betas_strat, indx, out_pre, type) {

  betas_df$BETA_Random <- betas_random
  betas_df$BETA_Strat <- betas_strat

  fwrite(betas_df, paste0(out_pre, "betas/genos-gwas_",indx, type), row.names=F,quote=F,sep="\t", col.names = F)
}

# Function to Calculate Va
calc_Va <- function(freq, betas_random, betas_strat, indx, va_name) {

  # Calulate Va = 2 * sum(B^2 * p * (1-p))
  Va_random <- 2 * sum((betas_random)^2 * freq * (1 - freq))
  Va_strat <- 2 * sum((betas_strat)^2 * freq * (1 - freq))

  # Create output df
  out <- as.data.frame(matrix(c(Va_random, Va_strat), nrow = 1))
  colnames(out) <- c("Va_Random", "Va_Strat")
  fwrite(out,paste0(out_pre, "Va/Va_",indx, va_name) ,row.names=F,quote=F,sep="\t", col.names = T)
}

# Function to generate "num" numbers of flipped effect sizes
new_betas <- function(beta_file, out_pre, type, freq, va_name) {

  betas_df <- fread(beta_file)
  colnames(betas_df) <- c("ID", "A1", "BETA_Random", "BETA_Strat")
  random_betas <- betas_df$BETA_Random
  strat_betas <- betas_df$BETA_Strat

  mat_random <- replicate(num, random_betas)
  mat_strat <- replicate(num, strat_betas)

  mat_random <- apply(mat_random, 2, flip)
  mat_strat <- apply(mat_strat, 2, flip)

  # Select only sites included in PRS
  allele_freq <- dplyr::inner_join(freq, betas_df, by = "ID")
  allele_freq <- allele_freq$ALT_FREQS

  for(i in 1:num) {
    write_file(betas_df, mat_random[,i], mat_strat[,i], i, out_pre, type)
    calc_Va(allele_freq, mat_random[,i], mat_strat[,i], i, va_name)
  }
}

# Load frequency file
freq <- fread(freq_file)

# General Resampled Betas for all types of PGS
new_betas(c_file, out_pre, ".c.betas", freq, ".c" )
new_betas(cp_file, out_pre, ".c.p.betas", freq, ".c.p")
new_betas(nc_file, out_pre, ".nc.betas", freq, ".nc")
new_betas(c_Tm_file, out_pre, ".c.betas-Tm", freq, ".c-Tm")
new_betas(cp_Tm_file, out_pre, ".c.p.betas-Tm", freq, ".c.p-Tm")
new_betas(nc_Tm_file, out_pre, ".nc.betas-Tm", freq, ".nc-Tm")









