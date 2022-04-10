## This script calculates Tm fromm plink output

## Requires: .eigenvec, .eigenval, .sscore, Tvec.txt

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_Tm.R <eigenvecs> <eigenvals> <sscore> <Tvec> <outfile Tm name> <outfile weights name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))

vecs_file = args[1] # eigenvectors
vals_file = args[2] # eigenvalues
proj_file = args[3] # .sscore (project of all test eigenvectors into eigenvalues)
tvec_file = args[4] # standardized test vector
id_file = args[5]
out_file = args[6] # Name for Tm file
out_file_weight=args[7] # Name for weights file

# Load test eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,2:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)


# Load test vector
std.tvec <- fread(tvec_file)
std.tvec$Tvec <- as.numeric(std.tvec$Tvec)


# Get the weights of each eigenvector
B <- t(vecs) %*% std.tvec$Tvec
Br = abs(B) / sum(abs(B))
out <- as_tibble(cbind(B, Br))
colnames(out) <- c("PC_weights", "Relative_PC_weights")
fwrite(out, out_file_weight,row.names=F,quote=F,sep="\t", col.names = T)

# Load Projected eigenvectors
proj_vecs <- fread(proj_file)
proj_vecs <- proj_vecs[,4:ncol(proj_vecs)]
proj_vecs <- apply(proj_vecs, 2, as.numeric)

# Load Eigenvalues
vals <- fread(vals_file)
vals <- vals$V1

# Calculate Tm
Tm <- proj_vecs %*% diag(1/sqrt(vals)) %*% B
Tm <- -2 * Tm

# Add to pop info
pops <- fread(id_file, header = T)
pops$Tm <- Tm

# Save Tm
fwrite(pops, out_file,row.names=F,quote=F,sep="\t", col.names = T)

