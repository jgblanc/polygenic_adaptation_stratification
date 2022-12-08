# This script takes in GWAS summary stats from plink2 and bins them into deciles based on effect size and averages effect size per decile

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript bin_avg_effect_sizes.R <uncorrect gwas> <tm gwas> <ID gwas> <out uncorrected> <out Tm> <out ID>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gwas_u = args[1] # Uncorrected summary stats
gwas_Tm = args[2] # Tm summary stats
gwas_ID = args[3] # ID summary stats
out_u = args[4] # Outpath for uncorrected ss
out_Tm = args[5] # Outpath for Tm ss
out_ID = args[6] # Outpath for ID ss

# Read in summary stats
gwas_u <- fread(gwas_u)
gwas_Tm <- fread(gwas_Tm)
gwas_ID <- fread(gwas_ID)

# Function to get the effect size for the T allele
flip_effect = function(gwas_df){
  gwas_df <- gwas_df %>% mutate(BETA = case_when(A1 == "A" ~ -BETA, A1== "T" ~ BETA))
  gwas_df$A1="T"
  return(gwas_df)
}

# Flip effect sizes to T allele
gwas_u <- flip_effect(gwas_u)
gwas_Tm <- flip_effect(gwas_Tm)
gwas_ID <- flip_effect(gwas_ID)

# Bin effect sizes into deciles and take average within each bin
bin_u <- gwas_u %>% mutate(es_bin = ntile(BETA, n=10)) %>% group_by(es_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))
bin_Tm <- gwas_Tm %>% mutate(es_bin = ntile(BETA, n=10)) %>% group_by(es_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))
bin_ID <- gwas_ID %>% mutate(es_bin = ntile(BETA, n=10)) %>% group_by(es_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))

# Save output
fwrite(bin_u, out_u, col.names=T,row.names=F,quote=F,sep="\t")
fwrite(bin_Tm, out_Tm, col.names=T,row.names=F,quote=F,sep="\t")
fwrite(bin_ID, out_ID, col.names=T,row.names=F,quote=F,sep="\t")

