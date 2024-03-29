# This Script uses the normal approximation to generate 2 sets of genotypes
library(data.table)

args=commandArgs(TRUE)
#rep <- args[1]
rep <- "A486"

# Set parameters
n <- 200
L_gwas <- 500
L <- 10000
ms <- 0

# Function to make Genotype matrix
make_genotype_matrix <- function(p, n, L) {
  G <- matrix(NA, nrow = n, ncol = L)
  for (i in 1:length(p)) {
    G[,i] <- rbinom(n,2,p[i])
  }
  return(G)
}


# Draw alleles from Beta
p_start <- rbeta(L, 1,3)
p_start_gwas <- rbeta(L_gwas, 1,3)

# Evolve GWAS alleles (ms = mean shift to add direction to frequency change)
p1GWAS <- p_start_gwas + rnorm(L_gwas, ms, sqrt((p_start_gwas)*(1-p_start_gwas)* 0.01))
p2GWAS <- p_start_gwas + rnorm(L_gwas, -ms, sqrt((p_start_gwas)*(1-p_start_gwas) *0.01))
rm <- which(p1GWAS < 0 | p1GWAS > 1 | p2GWAS < 0 | p2GWAS > 1)
if (length(rm) > 0) {
  p1GWAS <- p1GWAS[-rm]
  p2GWAS <- p2GWAS[-rm]
}

# Evolve neautral alleles vias drift
p1N <- p_start + rnorm(L, 0, sqrt((p_start_gwas)*(1-p_start_gwas)* 0.01))
p2N <- p_start + rnorm(L, 0, sqrt((p_start_gwas)*(1-p_start_gwas)* 0.01))
rm <- which(p1N < 0 | p1N > 1 | p2N < 0 | p2N > 1)
if (length(rm) > 0) {
  p1N <- p1N[-rm]
  p2N <- p2N[-rm]
}

# Make Genotype matrix for GWAS alleles by drawing from population
i1GWAS <- make_genotype_matrix(p1GWAS, n/2, length(p1GWAS))
i2GWAS <- make_genotype_matrix(p2GWAS, n/2, length(p2GWAS))

# Make Genotype matrix for neautral alleles by drawing from population
i1N <- make_genotype_matrix(p1N, n/2, length(p1N))
i2N <- make_genotype_matrix(p2N, n/2, length(p2N))

# Neutral Gentotype Matrix
X <- rbind(i1N, i2N)
if (length(which(colSums(X) == 0)) > 0 | length(which(colSums(X) == nrow(X))) > 0 ){
  if (length(which(colSums(X) == 0)) > 0) { X <- X[, -which(colSums(X) == 0)]}
  if (length(which(colSums(X) == nrow(X))) > 0) {X <- X[, -which(colSums(X) == nrow(X))]}
}

# GWAS Xentotype Matrix
X_gwas <- rbind(i1GWAS, i2GWAS)
if (length(which(colSums(X_gwas) == 0)) > 0 | length(which(colSums(X_gwas) == nrow(X_gwas))) > 0 ){
  if (length(which(colSums(X_gwas) == 0)) > 0) { X_gwas <- X_gwas[, -which(colSums(X_gwas) == 0)]}
  if (length(which(colSums(X_gwas) == nrow(X_gwas))) > 0) {X_gwas <- X_gwas[, -which(colSums(X_gwas) == nrow(X_gwas))]}
}

# Simulate Beta Hats
Bhat <- rnorm(ncol(X_gwas), 0,1)

# Compute Individual PGS
Z <- rowSums(X_gwas %*% diag(Bhat))

# Make test vector and compute Ztest
Tvec <- c(rep(1,(n/2)), rep(-1,(n/2)))

## Caclulate Eigenvec+Eigenval files

# Compute GRM
inv_het <- 1 / (2 * colMeans(X/2) * (1 - colMeans(X/2)))
Tmat <- matrix(-1/n, ncol = n, nrow = n)
diag(Tmat) <- (n-1)/n
cov_mat <- (1/(ncol(X)-1)) * Tmat %*% t(t(X) * inv_het) %*% t(X) %*% t(Tmat)

# Do eigen decompostion
myE <- eigen(cov_mat)
vecs <- myE$vectors
vals <- myE$values

# Format output
IDS <- seq(1,n)
vecs_df <- as.data.frame(cbind(IDS, IDS, vecs[,1:(n-1)]))
names <- rep("NA", n-1)
for (i in 1:(n-1)) {
  names[i] <- paste0("PC", i)
}
colnames(vecs_df) <- c("#FID", "IID", names)

# Save eigen
fwrite(vecs_df, paste0("../../output/Calculate_Tm/4PopSplit/", rep, "/C1/pca.eigenvec") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(vals[1:(n-1)]), paste0("../../output/Calculate_Tm/4PopSplit/", rep, "/C1/pca.eigenval") ,row.names=F,quote=F,sep="\t", col.names = F)

# Save test vector
fwrite(as.data.frame(Tvec), paste0("../../output/Calculate_Tm/4PopSplit/", rep, "/C1/Tvec.txt") ,row.names=F,quote=F,sep="\t", col.names = F)

# Save effect sizes
IDS <- seq(1, length(Bhat))
Allele <- rep("T", length(Bhat))
beta_df <- as.data.frame(cbind(IDS, Allele, Bhat, Bhat))
colnames(beta_df) <- c("ID", "A1", "BETA_Random", "BETA_strat")
fwrite(as.data.frame(beta_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-gwas_common.c.betas") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(beta_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-gwas_common.c.p.betas") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(beta_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-gwas_common.nc.betas") ,row.names=F,quote=F,sep="\t", col.names = T)

# Save allele frequency
freq <- colMeans(X_gwas/2)
freq_df <- as.data.frame(cbind(IDS, IDS, Allele, Allele, freq, rep(n, length(freq))))
colnames(freq_df) <- c("#CHROM", "ID", "REF", "ALT", "ALT_FREQS", "OBS_CT")
fwrite(as.data.frame(freq_df), paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-test_common.afreq"),row.names=F,quote=F,sep="\t", col.names = T)

# Save PGS
in_ID <- seq(1, n)
pgs_df <- as.data.frame(cbind(in_ID, rep(0,n), Z, Z))
colnames(pgs_df) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_SUM", "SCORE2_SUM")
fwrite(as.data.frame(pgs_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-test_common.c.sscore") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(pgs_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-test_common.c.p.sscore") ,row.names=F,quote=F,sep="\t", col.names = T)
fwrite(as.data.frame(pgs_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-test_common.nc.sscore") ,row.names=F,quote=F,sep="\t", col.names = T)

# Save "true" breeding valeu ("0")
in_ID <- seq(1, n)
true_df <- as.data.frame(cbind(in_ID, rep(0,n), rep(0,n)))
colnames(true_df) <- c("#IID", "NAMED_ALLELE_DOSAGE_SUM", "SCORE1_SUM")
fwrite(as.data.frame(true_df), paste0("../../output/PRS/4PopSplit/", rep, "/C1/h2-0/genos-test_common.true.sscore") ,row.names=F,quote=F,sep="\t", col.names = T)


