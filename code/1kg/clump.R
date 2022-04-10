args=commandArgs(TRUE)

if(length(args)<3){stop("<causal effects> <gwas results> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

ce_file = args[1]
gwas_file = args[2]
outfile = args[3]

# Read in causal sites
causal_effects <- fread(ce_file)
colnames(causal_effects) <- c("ID", "ALT", "CE")

# Read in GWAS results
gwas <- fread(gwas_file)
df <- inner_join(gwas, causal_effects)

# Function to get the effect size for the ALT allele
flip_effect = function(gwas_df,beta_colname){
  out <- gwas_df %>% mutate(BETA_Flip = case_when(ALT == A1 ~ BETA, ALT != A1 ~ -1 * BETA))
  return(out)
}

df_flipped <- flip_effect(df) %>% select("ID", "ALT", "BETA_Flip")
colnames(df_flipped)[3] <- "BETA"

fwrite(df_flipped,
       outfile,
       col.names=T,row.names=F,quote=F,sep="\t")
