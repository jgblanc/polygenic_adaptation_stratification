#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file> <output file prefix> <maximum PC>") }

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(R.utils)
}))

geffects_file=args[1] # Causal SNP IDs
gwas_file_prefix=args[2] # GWAS file prefix
output_file_prefix=args[3] # Out file prefix
max_pc = as.numeric(args[4])

# Function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  gwas_df = gwas_df[A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

# Function to select lowest p-value per chromosome
fclump <- function(df, CHR) {

  # Select only SNPS under threshold
  df <-  df %>% filter(CHROM == CHR)

  # Select lowest p-value
  min_p <- df %>% slice_min(P, with_ties = F)

  return(min_p)
}


## main function to get causal SNP effect sizes

main_casual <- function(gwas_path, out_path) {

  # Read in list of causal variants and clean data frame
  causal=fread(geffects_file)
  colnames(causal)=c("rsid","allele","esize")
  causal=causal%>%
    separate(rsid,into=c("CHROM","POS","ref","alt"),sep="_",remove=F)
  causal$POS=as.numeric(causal$POS)
  causal$CHROM=as.numeric(causal$CHROM)

  # Read in GWAS results
  gwas1=fread(paste(gwas_path,".pheno_strat",".glm.linear",sep=""),fill=T)
  colnames(gwas1)[1]="CHROM"

  # Flip to get correct effect sizes
  gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

  # Select effect sizes for causal variants
  gwas1.1 = gwas1[ID%in%causal$rsid]
  gwas.causal=gwas1.1[,c("ID","A1","BETA1")]
  colnames(gwas.causal) <- c("ID", "A1", "BETA_Strat")

  # Write Beta hat for causal effects
  fwrite(gwas.causal,
         paste(out_path,".c.betas",sep=""),
         col.names=T,row.names=F,quote=F,sep="\t")
}

#################################

main_nc <- function(gwas_path, out_path) {

  # Read in GWAS results
  gwas1=fread(paste(gwas_path,".pheno_strat",".glm.linear",sep=""),fill=T)
  colnames(gwas1)[1]="CHROM"

  # Flip to get correct effect sizes
  gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

  tmp = gwas1 %>% mutate(ID2 = ID) %>% separate(ID2, c("chr", "pos", "a1", "a2"), "_")
  gwas1= tmp %>% mutate(CHROM = chr) %>% select("CHROM", "POS", "ID", "A1", "BETA1", "P")
  nchrms = unique(gwas1$CHROM)
  gwas.red = fclump(gwas1, 1)
  gwas.red = gwas.red[order(gwas.red), ]

  # Repeat clumping for each chromosome
  for (i in 2:length(nchrms)) {
    new = fclump(gwas1, nchrms[i])
    new = new[order(new), ]
    gwas.red = rbind(gwas.red, new)
  }

  gwas.red = gwas.red %>% select(ID, A1, BETA1)
  gwas.red = gwas.red[,c("ID","A1","BETA1")]
  colnames(gwas.red) <- c("ID", "A1","BETA_Strat")


  # Write Beta hat from LD clumping
  fwrite(gwas.red,
         paste(out_path,".nc.betas",sep=""),
         col.names=T,row.names=F,quote=F,sep="\t")
}

#### Run for all types of GWASs

## Uncorrected
main_casual(gwas_file_prefix, output_file_prefix)
main_nc(gwas_file_prefix, output_file_prefix)

## Tm corrected
main_casual(paste0(gwas_file_prefix, "-Tm"), paste0(output_file_prefix, "-Tm"))
main_nc(paste0(gwas_file_prefix, "-Tm"), paste0(output_file_prefix, "-Tm"))

## ID corrected
main_casual(paste0(gwas_file_prefix, "-ID"), paste0(output_file_prefix, "-ID"))
main_nc(paste0(gwas_file_prefix, "-ID"), paste0(output_file_prefix, "-ID"))

## PC corrected
for (i in 1:max_pc) {

  main_casual(paste0(gwas_file_prefix, "-", i), paste0(output_file_prefix, "-", i))
  main_nc(paste0(gwas_file_prefix, "-", i), paste0(output_file_prefix, "-", i))

}




