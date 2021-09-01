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

# Load test vectors
tvecs <- fread(tvec_file)

# Get the weights of each eigenvector
B1 <- t(vecs) %*% tvecs$T1
B2 <- t(vecs) %*% tvecs$T2
B3 <- t(vecs) %*% tvecs$T3
B <- as.data.frame(cbind(B1,B2,B3))
colnames(B) <- c("T1", "T2", "T3")
fwrite(B, out_file_weight,row.names=F,quote=F,sep="\t", col.names = T)

# Load Projected eigenvectors
proj_vecs <- fread(proj_file)
proj_vecs <- proj_vecs[,5:ncol(proj_vecs)]
proj_vecs <- apply(proj_vecs, 2, as.numeric)

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Tm
Tm1 <- proj_vecs %*% diag(1/sqrt(vals)) %*% B1
Tm2 <- proj_vecs %*% diag(1/sqrt(vals)) %*% B2
Tm3 <- proj_vecs %*% diag(1/sqrt(vals)) %*% B3

# Save Tm
df <- as.data.frame(cbind(Tm1,Tm2,Tm3))
colnames(df) <- c("Tm1", "Tm2", "Tm3")
fwrite(df, out_file,row.names=F,quote=F,sep="\t", col.names = T)

