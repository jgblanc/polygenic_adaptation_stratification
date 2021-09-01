## This script calculates lambda_T

## Requires: .eigenvec .eigenval .tvec

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
tvecs <- fread(tvec_file)

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Lambda T
lambda_T1 <- t(tvecs$T1) %*% vecs %*% diag(vals) %*% t(vecs) %*% tvecs$T1
lambda_T2 <- t(tvecs$T2) %*% vecs %*% diag(vals) %*% t(vecs) %*% tvecs$T2
lambda_T3 <- t(tvecs$T3) %*% vecs %*% diag(vals) %*% t(vecs) %*% tvecs$T3

# Save Lambda T
lambda_T <- as.data.frame(c(lambda_T1, lambda_T2, lambda_T3))
fwrite(lambda_T, out_file,row.names=F,quote=F,sep="\t", col.names = T)
