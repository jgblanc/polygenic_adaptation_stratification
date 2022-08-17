# Split by CHR
# This script reads in the GWAS file and splits it into a separate file per chromosome so that the genotype matrix can be read into R
args=commandArgs(TRUE)
print(args)

if(length(args)!=4){stop("Rscript split_CHR.R <uncorrected GWAS>  <Tm GWAS> <ID GWAS> <out prefix>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(tidyverse)
}))

gwas_uncorrected = args[1]
gwas_Tm = args[2]
gwas_ID = args[3]
out_prefix = args[4]

# Function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  gwas_df = gwas_df[A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

# Read in files
gU <- fread(gwas_uncorrected)
colnames(gU)[1]="CHROM"
gTm <- fread(gwas_Tm)
colnames(gTm)[1]="CHROM"
gID <- fread(gwas_ID)
colnames(gID)[1]="CHROM"

# Flip effect sizes
gU <- flip_effect(gU,beta_colname = "BETA_strat")
gTm <- flip_effect(gTm,beta_colname = "BETA_strat")
gID <- flip_effect(gID,beta_colname = "BETA_strat")

# Separate ID file
gU_split <- gU %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_",remove =F)
gTm_split <- gTm %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_",remove=F)
gID_split <- gID %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_", remove=F)

# Loop through chromosomes and save files

# Uncorrected
maxChr <- max(as.numeric(gU_split$chr))
print(maxChr)
for (i in 0:(maxChr-1)) {
  df <- gU_split %>% filter(chr == (i+1)) %>% select(c("ID", "A1", "BETA_strat"))
  outfile <- paste0(out_prefix, i, ".betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}

# Tm
maxChr <- max(as.numeric(gTm_split$chr))
for (i in 0:(maxChr-1)) {
  df <- gTm_split %>% filter(chr == (i+1)) %>% select(c("ID", "A1", "BETA_strat"))
  outfile <- paste0(out_prefix, i, "-Tm.betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}

# ID
maxChr <- max(as.numeric(gID_split$chr))
for (i in 0:(maxChr-1)) {
  df <- gID_split %>% filter(chr == (i+1)) %>% select(c("ID", "A1", "BETA_strat"))
  outfile <- paste0(out_prefix, i, "-ID.betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}






