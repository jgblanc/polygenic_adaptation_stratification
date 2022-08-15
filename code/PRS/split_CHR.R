# Split by CHR
# This script reads in the GWAS file and splits it into a separate file per chromosome so that the genotype matrix can be read into R

if(length(args)!=4){stop("Rscript split_CHR.R <uncorrected GWAS>  <Tm GWAS> <ID GWAS> <out prefic>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gwas_uncorrected = args[1]
gwas_Tm = args[2]
gwas_ID = args[3]
out_prefix = args[4]

# Read in file
gU <- fread(gwas_uncorrected)
gTm <- fread(gwas_Tm)
gID <- fread(gwas_ID)

# Separate ID file
gU_split <- gU %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_")
gTm_split <- gTm %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_")
gID_split <- gID %>% separate(ID, c("chr", "tmp1", "tmp2", "tmp3"), "_")

# Loop through chromosomes and save files

# Uncorrected
maxChr <- max(as.numeric(gU_split$chr))
for (i in 1:maxChr) {
  df <- gU_split %>% filter(chr == i) %>% select(-c("tmp1","tmp2","tmp3","chr"))
  outfile <- paste0(out_prefix, i, ".betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}

# Tm
maxChr <- max(as.numeric(gTm_split$chr))
for (i in 1:maxChr) {
  df <- gTm_split %>% filter(chr == i) %>% select(-c("tmp1","tmp2","tmp3","chr"))
  outfile <- paste0(out_prefix, i, "-Tm.betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}

# ID
maxChr <- max(as.numeric(gID_split$chr))
for (i in 1:maxChr) {
  df <- gID_split %>% filter(chr == i) %>% select(-c("tmp1","tmp2","tmp3","chr"))
  outfile <- paste0(out_prefix, i, "-ID.betas")
  fwrite(df, outfile,row.names=F,quote=F,sep="\t", col.names = T)
}






