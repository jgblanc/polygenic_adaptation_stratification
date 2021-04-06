## This script calculates Tm fromm plink output

## Requires: .eigenvec, .eigenval, .sscore, Tvec.txt

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript calc_Tm.R <eigenvecs> <eigenvals> <tvec file> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

vecs_file = args[1] # eigenvectors
vals_file = args[2] # eigenvalues
tvec_file = args[3] # test vec file
out_file = args[4] # Name for Tm file

# Load test eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)

# Load test vector
std.tvec <- fread(tvec_file)
std.tvec <- std.tvec$V1

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Lambda T
lambda_T <- t(std.tvec) %*% vecs %*% diag(vals) %*% t(vecs) %*% std.tvec

# Save Lambda T
fwrite(lambda_T, out_file,row.names=F,quote=F,sep="\t", col.names = T)

