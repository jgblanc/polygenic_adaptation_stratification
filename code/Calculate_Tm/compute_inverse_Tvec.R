# The script uses the covariance matrix computed from plink (un-standardized allele counts) and computes (XX^T)^{-1}T

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript compute_inverse_Tvec.R <covariance matrix> <IDs> <Test Vector> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

cov_file = args[1]
ID_file = args[2]
Tvec_file = args[3]
outfile = args[4]


## Load in covariance matrix
cov_mat <- fread(cov_file)
cov_mat <- apply(cov_mat, 2, as.numeric)

## Do eigen decomposition on covariance matrix
myE <-eigen(cov_mat)

## Compute inverse covariance matrix
n <- ncol(cov_mat) - 1
FXXinv <- myE$vectors[,1:n] %*% diag(1/myE$values[1:n]) %*% t(myE$vectors[,1:n])

## Load Test vector
tvec <- fread(Tvec_file)
tvec <- tvec$Tvec

## Multiply inverse by Tvec
Tvec_inv <- FXXinv %*% tvec

## Load IDs and format output
IDs <- fread(ID_file)
out <- cbind(IDs, tvec, Tvec_inv)
colnames(out) <- c("#FID", "IID", "Tvec", "InTvec")

## Save output 
fwrite(out, outfile,col.names=T,row.names=F,quote=F,sep="\t")



