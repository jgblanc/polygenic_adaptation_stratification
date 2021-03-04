## Calculate Tm (small)

print("Running")

args=commandArgs(TRUE)

if(length(args)<4){stop("Rscript simphenotype_ge.R <test file> <gwas file> <popfile> <output_file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(snpStats)
}))

test_file = args[1] #genetic values
gwas_file = args[2] # pop file from msprime simulation
pop_file = args[3]
output_file = args[4] #name of output file

# Load genotype matricies
GWAS <- snpStats::read.plink(gwas_file)
TEST <- snpStats::read.plink(test_file)

# Convert to numeric matrices
X <- as(TEST$genotypes,"numeric")
X <- scale(X, scale = F)
M <- as(GWAS$genotypes, "numeric")
M <- scale(M, scale = F)

# Make Test vector
pops <- fread(pop_file, header = F)
pop <- dplyr::inner_join(pops, TEST$fam, by = c("V1"= "member")) %>% select("V2", "V3")
test_pops <- unique(pop$V3)
pop <- pop %>% mutate(tvec = case_when(V3 == test_pops[1] ~ 1, V3 == test_pops[2] ~ 0))
tvec <- pop$tvec
ctvec <- (tvec-mean(tvec))
Tvec <- ctvec/sqrt(sum(ctvec^2))


## get test covariance matrix and eigendecomposition
test.cov <- X %*% t(X) / (ncol(X)-1)
test.cov <- X %*% t(X) / (nrow(X))
X <- scale(X)
eig <- eigen(test.cov)
vecs <- eig$vectors
vals <- eig$values
n <- length(vals)
head(vals)

s <- svd(X)
u <- s$u
d <- s$d
v <- s$v
head(d^2 / (ncol(X) - 1))




# Calculate Tm
K = (M %*% t(X)) / (ncol(X) -1)
Tm = K %*% vecs[,1:(n-1)] %*% diag(1/vals[1:(n-1)]) %*% t(vecs[,1:(n-1)]) %*% Tvec

#ev = vecs[,1:(n-1)] %*% diag(1/vals[1:(n-1)]) %*% t(vecs[,1:(n-1)])
#sv = u[,1:(n-1)] %*% diag(1/d[1:(n-1)]) %*% t(u[,1:(n-1)])
#t = u[,1:(n-1)] %*% diag(1/d[1:(n-1)])

# Format datafile
tmp <- dplyr::inner_join(pops, GWAS$fam, by = c("V1"= "member")) %>% select("V2", "V3")
df <- as.data.frame(cbind(tmp$V2, tmp$V2, Tm))
colnames(df) <- c("FID","IID", "Tm")

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)

