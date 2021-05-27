# This script takes actual estimated effect sizes and randomly flips the effect sizes to generate a given number of resampled effect sizes to general an emprical null


args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_Tm.R <c.betas> <c.p.betas> <n.c.betas> <num resample> <output prefix>
                         <true.sscore> <Tvec.txt> <outfile name> <file with number of snps>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

c_file = args[1] # causal betas
#c_file = "output/PRS/4PopSplit/E1/C1/h2-0/env-0.0/genos-gwas_common.c.betas"
cp_file = args[2] # causal p-value betas
#cp_file = "output/PRS/4PopSplit/E1/C1/h2-0/env-0.0/genos-gwas_common.c.p.betas"
nc_file = args[3] # clumped betas
#nc_file = "output/PRS/4PopSplit/E1/C1/h2-0/env-0.0/genos-gwas_common.nc.betas"
num = args[4] # number of times to resample
out_pre = args[5] # output prefix


# Function to flip effect sizes
flip <- function(betas, num) {
  new_betas <- sample(c(-1,1), num,  replace = T) * betas
  return(new_betas)
}

# Function to write new betas to a file
write_file <- function(betas_random, betas_strat, indx, out_pre, type) {

  betas_df <- fread(beta_file)
  print(betas_random)
  betas_df$BETA_Random <- betas_random
  betas_df$BETA_Strat <- betas_strat

  fwrite(betas_df, paste0(out_pre, indx, type))
}

# Function to generate "num" numbers of flipped effect sizes
new_betas <- function(beta_file, out_pre, type) {

  betas_df <- fread(beta_file)
  random_betas <- betas_df$BETA_Random
  strat_betas <- betas_df$BETA_Strat

  mat_random <- replicate(num, random_betas)
  mat_strat <- replicate(num, strat_betas)

  mat_random <- apply(mat_random, 2, flip, num=num)
  mat_strat <- apply(mat_strat, 2, flip, num=num)

  for(i in 1:num) {
    write_file(mat_random[,i], mat_strat[,i], i, out_pre)
  }
}

new_betas(c_file, out_pre, ".c.betas")
new_betas(cp_file, out_pre, ".c.p.betas")
new_betas(nc_file, out_pre, ".nc.betas")









