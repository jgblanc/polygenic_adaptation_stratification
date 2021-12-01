## This script calculates Tm fromm plink output

## Requires: .eigenvec, .eigenval, .sscore, Tvec.txt

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript calc_Tm.R <eigenvecs> <eigenvals> <sscore> <Tvec> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

vecs_file = args[1] # eigenvectors
vals_file = args[2] # eigenvalues
proj_file = args[3] # .sscore (project of all test eigenvectors into eigenvalues)
tvec_file = args[4] # standardized test vector
out_file = args[5] # Name for Tm file
out_file_weight=args[6] # Name for wieghts file

# Load test eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)
#head(vecs)

# Load test vector
std.tvec <- fread(tvec_file)
std.tvec$V1 <- as.numeric(std.tvec$V1)
#print(std.tvec)
#print(str(std.tvec))

# Get the weights of each eigenvector
B <- t(vecs) %*% std.tvec$V1
#fwrite(as.data.frame(B[,1]), out_file_weight,row.names=F,quote=F,sep="\t", col.names = T)

# Load Projected eigenvectors
proj_vecs <- fread(proj_file)
proj_vecs <- proj_vecs[,5:ncol(proj_vecs)]
proj_vecs <- apply(proj_vecs, 2, as.numeric)
#head(proj_vecs)

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Tm
Tm <- proj_vecs %*% diag(1/sqrt(vals)) %*% B
Tm <- -2 * Tm

# Save Tm
fwrite(Tm, out_file,row.names=F,quote=F,sep="\t", col.names = T)

