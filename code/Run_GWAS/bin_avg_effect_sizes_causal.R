# This script takes in GWAS summary stats from plink2 and bins them into deciles based on effect size and averages effect size per decile

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript bin_avg_effect_sizes.R <r> <uncorrect gwas> <tm gwas> <ID gwas> <out uncorrected> <out Tm> <out ID>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

r_file =args[1]
gwas_u = args[2] # Uncorrected summary stats
gwas_Tm = args[3] # Tm summary stats
gwas_ID = args[4] # ID summary stats
true_es = args[5]
out_u = args[6] # Outpath for uncorrected ss
out_Tm = args[7] # Outpath for Tm ss
out_ID = args[8] # Outpath for ID ss
out_te = args[9]

# Read in r bins
r_bin <- fread(r_file)

# Read in summary stats
gwas_u <- fread(gwas_u)
gwas_Tm <- fread(gwas_Tm)
gwas_ID <- fread(gwas_ID)

# Read in true effect sizes
true_es <- fread(true_es)
colnames(true_es) <- c("ID", "A1", "BETA")
ce_ids <- true_es$ID

# Filter to only include causal sites
gwas_u <- gwas_u %>% filter(ID %in% ce_ids)
gwas_Tm <- gwas_Tm %>% filter(ID %in% ce_ids)
gwas_ID <- gwas_ID %>% filter(ID %in% ce_ids)

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
print(nrow(gwas_u))

# Bin effect sizes into deciles and take average within each bin
bin_u <- inner_join(gwas_u, r_bin) %>% select(ID, r_bin, BETA) %>% group_by(r_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))
bin_Tm <- inner_join(gwas_Tm, r_bin) %>% select(ID, r_bin, BETA) %>% group_by(r_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))
bin_ID <- inner_join(gwas_ID, r_bin) %>% select(ID, r_bin, BETA) %>% group_by(r_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))
bin_te <- inner_join(true_es, r_bin) %>% select(ID, r_bin, BETA) %>% group_by(r_bin) %>% summarise(avg = mean(BETA), SD = sd(BETA))

# Save output
fwrite(bin_u, out_u, col.names=T,row.names=F,quote=F,sep="\t")
fwrite(bin_Tm, out_Tm, col.names=T,row.names=F,quote=F,sep="\t")
fwrite(bin_ID, out_ID, col.names=T,row.names=F,quote=F,sep="\t")
fwrite(bin_te, out_te, col.names=T,row.names=F,quote=F,sep="\t")

