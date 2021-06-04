#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <p-value threshold> <output file prefix>") }

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(R.utils)
}))

geffects_file=args[1]
#geffects_file="~/polygenic_adaptation_stratification/output/Simulate_Phenotypes/4PopSplit/S1/C2/h2-0/genos-gwas_common.effects.txt"
gwas_file_prefix=args[2]
#gwas_file_prefix="~/polygenic_adaptation_stratification/output/Run_GWAS/4PopSplit/S1/C2/h2-0/env-0.0/genos-gwas_common-Tm"
pval_threshold=as.numeric(args[3])
output_file_prefix=args[4]
#output_file_prefix="~/polygenic_adaptation_stratification/output/PRS/4PopSplit/V1/C1/genos-gwas_common"


#read in list of causal variants (true simulated effects)
causal=fread(geffects_file)
colnames(causal)=c("rsid","allele","esize")

causal=causal%>%
  separate(rsid,into=c("CHROM","POS","ref","alt"),sep="_",remove=F)

causal$POS=as.numeric(causal$POS)
causal$CHROM=as.numeric(causal$CHROM)

print("reading gwas files")
#gwas effects
gwas1=fread(paste(gwas_file_prefix,".pheno_strat",".glm.linear",sep=""),fill=T)

colnames(gwas1)[1]="CHROM"

#function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  #gwas_df = gwas_df[ ID %in% causal$rsid, .(CHROM,POS,ID,A1,BETA,P)]
  gwas_df = gwas_df[ A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[ A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}

gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

print("filtering causal variants")
#select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]



gwas.causal=gwas1.1[,c("ID","A1","BETA1")]
colnames(gwas.causal) <- c("ID", "A1", "BETA_Strat")

fwrite(gwas.causal,
       paste(output_file_prefix,".c.betas",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

print("filtering variants under a pvalue threshold")

#write function to select causal variants below some p-value threshold
fcausal_p = function(gwas_df,beta_colname,pvalue=pval_threshold){

  gwas_df[ P > pvalue, (beta_colname) := 0]
  return(gwas_df)
}

fcausal_p = function(df,pvalue=pval_threshold){

  df=df%>%
    filter(P < pvalue)
  return(df)
}

gwas1.2 = fcausal_p( gwas1.1,pval_threshold )
gwas.causal.p = gwas1.2[,c("ID","A1","BETA1")]

colnames(gwas.causal.p) <- c("ID", "A1", "BETA_Strat")

# Add dummy row with effect size zero
if (nrow(gwas.causal.p) < 1) {
  gwas.causal.p <- rbind(gwas.causal.p, gwas.causal[1,])
  gwas.causal.p[1,3] <- 0
  gwas.causal.p[1,4] <- 0
}

fwrite(gwas.causal.p,
       paste(output_file_prefix, ".c.p.betas" , sep=""),
       col.names=T, row.names=F , quote=F , sep="\t")

print("ld clumping")

fclump=function(df,pcutoff){
  if(missing(pcutoff)){
    df.red=df
  }else{
    df.red=df%>%
      filter(P<pcutoff)
  }

  df.red$window=NA
  for(i in 1:101){
    #start=((i-1)*1e5 - 5e4) +1
    start=((i-1)*1e5) +1
    stop=start+1e5
    df.red$window[which((df.red$POS>=start) & (df.red$POS<stop))]=i
  }

  df.red$window_name=paste(df.red$CHROM,df.red$window,sep="_")

  return(df.red)

}

flead=function(df){
  df.red=df%>%
    slice_min(P, with_ties = F)
  return(df.red)
}

gwas1.red=fclump(gwas1,pval_threshold)%>%
  group_by(window_name)%>%
  flead(.)%>%
  ungroup()%>%
  select(ID,A1,BETA1)


gwas.red = gwas1.red[,c("ID","A1","BETA1")]
colnames(gwas.red) <- c("ID", "A1","BETA_Strat")

# Add dummy row with effect size zero
if (nrow(gwas.red) < 1) {
  gwas.red <- rbind(gwas.red, gwas.causal[1,])
  gwas.red[1,3] <- 0
  gwas.red[1,4] <- 0
}

fwrite(gwas.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")
