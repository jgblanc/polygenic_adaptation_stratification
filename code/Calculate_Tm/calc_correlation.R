# This script calculates the correlation between test vector and Tm

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript calc_correlation.R <Tm> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

Tm_file = args[1] # Path to Tm file
#Tm_file = "output/Calculate_TGWAS/4PopSplit/S1/C1/Tm.txt"
outfile = args[2] # Path to output directory


# Read in Tm
Tm <- fread(Tm_file)

# Make fake test vector with block structure
Tvec <- as.data.frame(c(rep(0.5, nrow(Tm)/2), rep(-0.5, nrow(Tm)/2)))
colnames(Tvec) <- "Tvec"

# Compute correlation
out <- cor(Tm$Tm, Tvec$Tvec)
out <- as.data.frame(out)
colnames(out) <- "correlation"

# Write to file
fwrite(out, outfile,row.names=F,quote=F,sep="\t", col.names = T)
