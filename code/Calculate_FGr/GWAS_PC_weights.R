## This script calculates the relative weights of GWAS PCs on T^GWAS

## Requires: gwas_pca.eigenvec, Tm.txt

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_Tm.R <eigenvecs> <Tm> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))


vecs_file = args[1] # eigenvectors
Tm_file = args[2] # Tm
out_file = args[3]

# Load gwas eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)
vecs <- scale(vecs)


# Load T^GWAS
Tm = fread(Tm_file)
print(head(Tm))
Tm = Tm$Tm

# Normalize T to have var 1
Tm <- scale(Tm)

# Get weights
B = t(vecs) %*% Tm

# Get relative weights
Br = abs(B) / sum(abs(B))

# Convert to dataframe
B <- as_tibble(cbind(B, Br))
colnames(B) <- c("PC_weights", "Relative_PC_weights")

# Save weights
fwrite(B, out_file,row.names=F,quote=F,sep="\t", col.names = T)
