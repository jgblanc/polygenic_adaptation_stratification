# This script takes a setof effect sizes and test vector and computes Q

args=commandArgs(TRUE)

if(length(args)!=4){stop("Rscript calcQ_4PopSplit.R <prefix to effect sizes> <test panel genotypes prefix> <test vector>
                         <number of times to resample in empirical null> <output prefix> <true effect sizes> <list of PCs>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

gwas_prefix = args[1]
print(gwas_prefix)
test_prefix = args[2]
print(test_prefix)
out_file = args[3]
evec_file = args[4]


# Function to read in genotype matrix for a set of variants
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Read in genotype matrices
vars <- fread(paste0(gwas_prefix, ".pvar"), skip = 200)
G <- read_genos(gwas_prefix, vars)
X <- read_genos(test_prefix, vars)
print(dim(G))
print(dim(X))

# Read in PCs
PCs <- fread(evec_file)
u1 <- PCs$PC1

# Set parameters
prob <- 1
m <- nrow(G)
n <- nrow(X)
L <- ncol(G)
ncausal_list <- c(1, 2, 10, 100, 500, 1000, 2000)

#### Only once #####

# Get allele frequency diff in TP
pC <- colMeans(X[1:1000,])/2
pD <- colMeans(X[1001:2000,])/2
diff <- pD - pC

# Mean center both genotype matrices
G <- scale(G, scale = F)
X <- scale(X, scale = F)

# Make 2 pop test vector
tvec <- scale(c(rep(0,n/2),rep(1,n/2)))

# Compute FGr
r.all <- t(X) %*% tvec
Gvar <- apply(G, 2, var)
FGr <- G %*% (r.all/Gvar)
FGr <- scale(FGr)

########################

out <- matrix(NA, nrow = length(ncausal_list), ncol = 7)

for (j in 1:length(ncausal_list)) {

  print(j)

  ncausal <- ncausal_list[j]

  # Generate effect sizes so there is a true signal
  indx <- sample(seq(1, L), ncausal)
  diff.causal <- diff[indx]
  B <- rnorm(ncausal, 0, L/ncausal)
  for (i in 1:ncausal){
    if (diff.causal[i] >= 0) {
      B[i] <- sample(c(-1, 1),1, prob = c((1-prob), prob)) * abs(B[i])
    } else {
      B[i] <- sample(c(1, -1),1, prob = c((1-prob), prob)) * abs(B[i])
    }
  }
  betas <- rep(0, L)
  betas[indx] <- B

  # Calculate True GV
  gwas.gvs <- G[, indx]%*%as.matrix(B)
  test.gvs <- X[, indx]%*%as.matrix(B)

  # Simulate Phenotypes
  phenos <- gwas.gvs+rnorm(m)

  # Compute Bhat including Tm
  Tm_Bhat <- numeric()
  for(k in 1:L){
    Tm_Bhat[k] <- lm(phenos~G[,k] + FGr)$coef[2]
  }

  # Compute Bhat including PC1
  PC_Bhat <- numeric()
  for(k in 1:L){
    PC_Bhat[k] <- lm(phenos~G[,k] + u1)$coef[2]
  }

  # Compute different versions of q
  r.causal <- r.all[indx]
  true.q <- (1/(n-1)) * (B %*% r.causal)
  S.q <- (1/(n-1)) * (Tm_Bhat[indx] %*% r.causal)
  L.q <- (1/(n-1)) * (Tm_Bhat %*% r.all)
  eq <- true.q * (1 - (ncausal/L))
  PC_S.q <- (1/(n-1)) * (PC_Bhat[indx] %*% r.causal)
  PC_L.q <- (1/(n-1)) * (PC_Bhat %*% r.all)

  tmp <- c(true.q, S.q, L.q, ncausal, eq, PC_S.q, PC_L.q)
  out[j,] <- tmp
}

# Save output
print(out)
fwrite(out, out_file,row.names=F,quote=F,sep="\t", col.names = T)


