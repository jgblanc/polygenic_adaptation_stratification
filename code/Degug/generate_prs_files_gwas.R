# This Script uses msprime genotype matrices to calculate Qx
library(data.table)
library(pgenlibr)
library(dplyr)
library(stringr)

setwd("~/polygenic_adaptation_stratification/code/Degug/")
args=commandArgs(TRUE)
rep <- args[1]
rep <- "M1"

seed <- str_split(rep, "")
set.seed(as.numeric(seed[[1]][2]))

# Set parameters
L_gwas <- 500

num_f <- fread(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-gwas_common.afreq"))
num_var <- nrow(num_f)

# Read in test panel
pvar <- NewPvar(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-test_common.pvar"))
pvar <- NewPvar(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-test_common.pvar"))
d1 <- NewPgen(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-test_common.pgen"))
X <- ReadList(d1,seq(1,as.numeric(num_var)), meanimpute=F)
n <- nrow(X)
L <- ncol(X)

# Read in gwas panel
pvar <- NewPvar(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-gwas_common.pvar"))
d1 <- NewPgen(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-gwas_common.pgen"))
X_gwas <- ReadList(d1,seq(1,as.numeric(num_var)), meanimpute=F)

# Read in GWAS effect sizes
betas <- fread(paste0("../../output/Run_GWAS/4PopSplit/", rep, "/C1/h2-0/env-0/genos-gwas_common.pheno_random.glm.linear"))

# Sample L_gwas number effect sizes
gwas_indx <- sample(seq(0, ncol(X_gwas)), L_gwas)
X_gwas <- X_gwas[, gwas_indx]
Bhat <- betas[gwas_indx, 9]
Bhat <- Bhat$BETA

# Compute Individual PGS
Z <- rowSums(X_gwas %*% diag(Bhat))

# Make test vector

# Read in Fam file
fam <- fread(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/C1/genos-gwas_common.psam"))
colnames(fam) <- c("IID", "FID", "SEX")

# Make standardized test vec
pops <- fread(paste0("../../output/Simulate_Genotypes/4PopSplit/", rep, "/genos.pop"), header = F)
pop <- dplyr::inner_join(pops, fam, by = c("V1"= "IID")) %>% select("V2", "V3")
test_pops <- unique(pop$V3)
pop <- pop %>% mutate(tvec = case_when(V3 == test_pops[1] ~ 1, V3 == test_pops[2] ~ -1))
Tvec <- pop$tvec

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
