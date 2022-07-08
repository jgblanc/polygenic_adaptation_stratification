# This script calculates the correlation between test vector and Tm

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_correlation.R <Tm> <Tvec> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_file = args[1] # Path to Tm file
#Tm_file = "output/Calculate_TGWAS/4PopSplit/S1/C1/Tm.txt"
tvec_file = args[2] # Path to test vectors
#tvec_file = "output/Calculate_TGWAS/4PopSplit/S1/C1/Tvec.txt"
outfile = args[3] # Path to output directory


# Read in Tm
Tm <- fread(Tm_file)

# Read in Tvec and downsample to same size as Tm
Tvec <- fread(tvec_file)
Tvec <- Tvec %>% group_by(Tvec) %>% sample_n(nrow(Tm)/2)

# Compute correlation
out <- cor(Tm$Tm, Tvec$Tvec)
out <- as.data.frame(out)
colnames(out) <- "correlation"

# Write to file
fwrite(out, outfile,row.names=F,quote=F,sep="\t", col.names = T)
