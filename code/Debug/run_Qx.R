library(data.table)
library(dplyr)
library(ggplot2)

# Note these are small msprime simulations with 200 samples per panel and only 2 chromosomes

rep <- "M3"

# Read in GWAS Panel Genotypes
X_gwas <- readRDS("~/Desktop/GWAS.rds")
n <- nrow(X_gwas) # Number of Individuals
L <- ncol(X_gwas) # Number of total loci (no filtering except at least one copy per panel)

# Draw Random Phenotype for each individual (this matches the phenotypes fed into pipeline)
set.seed(30)
pheno <- scale(rnorm(n,0,1), scale = T)
df <- as.data.frame(cbind(pheno, Tvec <- c(rep(-1, 100)/100,rep(1, 100)/100) * (1/2)))
ggplot(df, aes(x=V1, fill = as.character(V2))) + geom_histogram(position="identity")

# Do GWAS by repgressing phenotype on genotype (these match plink output once you flip to the right allele) (take like a minute to run)
est_betas <- rep(0, L)
for(i in 1:L){
  mod <- lm(pheno ~ X_gwas[,i])
  est_betas[i] <- mod$coefficients[2]
}

# Pick 500 random sites to include in the PRS and get their estimated effect sizes
set.seed(3)
gwas_indx <- sample(seq(0, ncol(X_gwas)), 500)
Bhat <- est_betas[gwas_indx]
Bhat <- rnorm(500, 0, 0.05) #Uncomment this line to re-run with randomly drawn Beta hats


####################################

# Read in Test Panel Genotypes
X <- readRDS("~/Desktop/TEST.rds")

# Compute Individual PGS
Z <- rowSums(X[,gwas_indx] %*% diag(Bhat))
#Z <- X[, gwas_indx] %*% Bhat

# Make Test Vector (in this example I know that there are exactly 100 indviduals in pop B and 100 in pop D)
Tvec <- c(rep(-1, 100)/100,rep(1, 100)/100) * (1/2)

# Compute Numerator of Qx (this is mean(PGS_B) - mean(PGS_D)/2)
Ztest <- Z %*% Tvec
#Ztest2 <- (mean(Z[1:100]) - mean(Z[101:200])) / 2

## Caclulate Eigenvec+Eigenval files
# Compute GRM
inv_het <- 1 / (2 * colMeans(X/2) * (1 - colMeans(X/2)))
Tmat <- matrix(-1/n, ncol = n, nrow = n)
diag(Tmat) <- (n-1)/n
cov_mat <- (1/(ncol(X))) * Tmat %*% t(t(X) * inv_het) %*% t(X) %*% t(Tmat)

# Do eigen decompostion - these match the plink output (I'm decomposing GRM and then rebuilding it to mimic the pipeline where we only have the eigenvec/vals and not the whole GRM)
myE <- eigen(cov_mat)
vecs <- myE$vectors
vals <- myE$values

K <- myE$vectors[,1:199] %*% diag(myE$values[1:199]) %*% t(myE$vectors[,1:199]) # Recontruct GRM
Fmat <-  (Tvec %*% K %*% Tvec) # This is 2*Fst

# Compute Va
Va <- 2*sum(Bhat^2 * (colMeans(X[,gwas_indx]/2) * (1 - colMeans(X[,gwas_indx]/2))))


# Compute Qx
Qx <- (t(Ztest) %*% solve(Fmat) %*% Ztest) / (Va)
print(Qx)

# When using random betas Qx=0.84
# When using gwas estimated betas Qx=54

