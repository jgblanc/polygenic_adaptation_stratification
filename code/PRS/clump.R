
#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <causal effects file> <prefix for glm.linear file (w/o phenotype or correction)> <p-value threshold> <output file prefix>") }

library(data.table)
library(dplyr)
library(tidyr)
library(R.utils)

geffects_file=args[1]
#geffects_file="~/polygenic_adaptation_stratification/output/Simulate_Phenotypes/4PopSplit/V1/C1/genos-gwas_common.effects.txt"
gwas_file_prefix=args[2]
#gwas_file_prefix="~/polygenic_adaptation_stratification/output/Run_GWAS/4PopSplit/V1/C1/genos-gwas_common"
pval_threshold=args[3]
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
gwas1=fread(paste(gwas_file_prefix,".pheno_random",".glm.linear",sep=""),fill=T)
gwas2=fread(paste(gwas_file_prefix, ".pheno_strat",".glm.linear",sep=""),fill=T)

colnames(gwas1)[1]=colnames(gwas2)[1]="CHROM"

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
gwas2 = flip_effect(gwas2,beta_colname = "BETA2")

print("filtering causal variants")
#select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]
gwas2.1 = gwas2[ID%in%causal$rsid]


gwas.causal=cbind( gwas1.1[,c("ID","A1","BETA1")],
                    gwas2.1[,c("BETA2")])
colnames(gwas.causal) <- c("ID", "A1", "BETA_Random", "BETA_Strat")

fwrite(gwas.causal,
       paste(output_file_prefix,".c.betas",sep=""),
       col.names=T,row.names=F,quote=F,sep="\t")

print("filtering variants under a pvalue threshold")

#write function to select causal variants below some p-value threshold
fcausal_p = function(gwas_df,beta_colname,pvalue=pval_threshold){

    gwas_df[ P > pvalue, (beta_colname) := 0]
    return(gwas_df)
}

gwas1.2 = fcausal_p( gwas1.1, "BETA1" )
gwas2.2 = fcausal_p( gwas2.1, "BETA2" )

gwas.causal.p = cbind( gwas1.2[,c("ID","A1","BETA1")],
                    gwas2.2[,c("BETA2")])

colnames(gwas.causal.p) <- c("ID", "A1", "BETA_Random", "BETA_Strat")

fwrite(gwas.causal.p,
       paste(output_file_prefix, ".c.p.betas" , sep=""),
       col.names=T, row.names=F , quote=F , sep="\t")

print("ld clumping")
#write function to assign each SNP to window, then find a single hit within each
fclump=function(gwas_df,pcutoff=pval_threshold){
 if(missing(pcutoff)){
   df.red=gwas_df
 }else{
   df.red = gwas_df[P<pcutoff]
 }


 for(i in 1:101){
   if(i == 1){
     start=0
     stop=start+1e5
   }else{

   start=((i-1)*1e5) +1
   stop=start + 1e5 -1
   }

   df.red[ (POS>=start &POS<stop) , window:=i]
 }
 df.red[,window_name:= paste(CHROM,window,sep="_") ]
 return(df.red)

}

flead=function(df){
  df.red = df[P == min(P)]
 return(df.red)
}

#causal=fclump(causal)
gwas1.red = fclump(gwas1,pval_threshold)
gwas1.red = gwas1.red[,flead(.SD),by=.(window_name)]
gwas1.red = unique(gwas1.red, by ="window_name")
gwas1.red = gwas1.red[,.(ID,A1,BETA1)]

gwas2.red = fclump(gwas2,pval_threshold)
gwas2.red = gwas2.red[,flead(.SD),by=.(window_name)]
gwas2.red = unique(gwas2.red, by ="window_name")
gwas2.red = gwas2.red[,.(ID,A1,BETA2)]


gwas.red = merge(gwas1.red[,c("ID","A1","BETA1")], gwas2.red[,c("ID","BETA2")], by="ID", all=TRUE)
gwas.red$A1="T"

gwas.red[is.na(BETA1),BETA1:=0]
gwas.red[is.na(BETA2),BETA2:=0]
colnames(gwas.red) <- c("ID", "A1", "BETA_Random", "BETA_Strat")

colnames(gwas.causal.p) <- c("ID", "A1", "BETA_Random", "BETA_Strat")
fwrite(gwas.red,
       paste(output_file_prefix,".nc.betas",sep=""),
       col.names=F,row.names=F,quote=F,sep="\t")
