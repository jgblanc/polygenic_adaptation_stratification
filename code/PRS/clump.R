#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file> <p-value threshold> <output file prefix>") }

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(R.utils)
}))

geffects_file=args[1] # Causal SNP IDs
gwas_file_prefix=args[2] # GWAS file prefix
pval_threshold=as.numeric(args[3]) # P-value threshold
output_file_prefix=args[4] # Out file prefix

# Read in list of causal variants and clean data frame
causal=fread(geffects_file)
colnames(causal)=c("rsid","allele","esize")
causal=causal%>%
  separate(rsid,into=c("CHROM","POS","ref","alt"),sep="_",remove=F)
causal$POS=as.numeric(causal$POS)
causal$CHROM=as.numeric(causal$CHROM)

# Read in GWAS results
gwas1=fread(paste(gwas_file_prefix,".pheno_strat",".glm.linear",sep=""),fill=T)
colnames(gwas1)[1]="CHROM"

# Function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  gwas_df = gwas_df[A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

# Flip to get correct effect sizes
gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

# Select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]
gwas.causal=gwas1.1[,c("ID","A1","BETA1")]
colnames(gwas.causal) <- c("ID", "A1", "BETA_Strat")

# Write Beta hat for causal effects
fwrite(gwas.causal,
       paste(output_file_prefix,".c.betas",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

##################################

# Write function to select causal variants below some p-value threshold
fcausal_p = function(df,pvalue=pval_threshold){

  df=df%>%
    filter(P < pvalue)
  return(df)
}

gwas1.2 = fcausal_p( gwas1.1,pval_threshold )
gwas.causal.p = gwas1.2[,c("ID","A1","BETA1")]
colnames(gwas.causal.p) <- c("ID", "A1", "BETA_Strat")

# Add dummy row with effect size zero if there are no variants under the threshold
if (nrow(gwas.causal.p) < 1) {
  gwas.causal.p <- rbind(gwas.causal.p, gwas.causal[1,])
  gwas.causal.p[1,3] <- 0
  gwas.causal.p[1,4] <- 0
}

# Write Beta hat for causal effects under a p-value threshold
fwrite(gwas.causal.p,
       paste(output_file_prefix, ".c.p.betas" , sep=""),
       col.names=T, row.names=F , quote=F , sep="\t")

#################################

# Function to select lowest p-value per chromosome
fclump <- function(df, pt, CHR) {

  # Select only SNPS under threshold
  df <-  gwas1 %>% filter(P < pt) %>% filter(CHROM == CHR)

  # Select lowest p-value
  min_p <- df %>% slice_min(P, with_ties = F)

  return(min_p)
}

tmp = gwas1 %>% mutate(ID2 = ID) %>% separate(ID2, c("chr", "pos", "a1", "a2"), "_")
gwas1= tmp %>% mutate(CHROM = chr) %>% select("CHROM", "POS", "ID", "A1", "BETA1", "P")
nchrms = unique(gwas1$CHROM)
gwas.red = fclump(gwas1, pval_threshold, 1)
gwas.red = gwas.red[order(gwas.red), ]

# Repeat clumping for each chromosome
for (i in 2:length(nchrms)) {
  new = fclump(gwas1, pval_threshold, nchrms[i])
  new = new[order(new), ]
  gwas.red = rbind(gwas.red, new)
}

gwas.red = gwas.red %>% select(ID, A1, BETA1)
gwas.red = gwas.red[,c("ID","A1","BETA1")]
colnames(gwas.red) <- c("ID", "A1","BETA_Strat")

# Add dummy row with effect size zero if there are no variants under the p-value threshold
if (nrow(gwas.red) < 1) {
  gwas.red <- rbind(gwas.red, gwas.causal[1,])
  gwas.red[1,3] <- 0
  gwas.red[1,4] <- 0
}

# Write Beta hat from LD clumping
fwrite(gwas.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
