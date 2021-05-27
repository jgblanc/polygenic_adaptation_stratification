## This script calculates Va for 3 types of GWASs

## Requires: test panel allele frequency, betas (c, c.p, and n.c)

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_Tm.R <freq> <c.betas> <c.p.betas> <n.c.betas> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

freq_file = args[1] # test panel freq
c_file = args[2] # causal beta estimates
cp_file = args[3] # causal beta under threshold estimates
nc_file = args[4] # clumped beta estimates
out_file = args[5] # outfile prefix

# Load frequency file
freq <- fread(freq_file)

# Load beta files
c_betas <- fread(c_file)
cp_betas <- fread(cp_file)
nc_betas <- fread(nc_file)

# Function to Calculate Va
calc_Va <- function(freq, betas) {

  # Select only sites included in PRS
  colnames(betas) <- c("ID", "ALT", "Beta_Random", "Beta_Strat")
  df <- dplyr::inner_join(freq, betas)

  # Calulate Va = 2 * sum(B^2 * p * (1-p))
  Va_random <- 2 * sum((df$Beta_Random)^2 * df$ALT_FREQS * (1 - df$ALT_FREQS))
  Va_strat <- 2 * sum((df$Beta_Strat)^2 * df$ALT_FREQS * (1 - df$ALT_FREQS))

  # Create output df
  out <- as.data.frame(matrix(c(Va_random, Va_strat), nrow = 1))
  colnames(out) <- c("Va_Random", "Va_Strat")
  return(out)
}

out_c <- calc_Va(freq, c_betas)
out_cp <- calc_Va(freq, cp_betas)
out_nc <- calc_Va(freq, nc_betas)

out <- rbind(out_c, out_cp, out_nc)
row.names(out) <- c("c", "c.p", "nc")

# Save Va's
fwrite(out, out_file,row.names=T,quote=F,sep="\t", col.names = T)
