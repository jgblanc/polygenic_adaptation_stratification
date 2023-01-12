# This script takes in selected SNPs, re-estimates effect sizes jointly, and computes q for different numbers of SNPs

args=commandArgs(TRUE)

if(length(args)<11){stop("Rscript joint_effect_sizes.R <asecertained snps> <causal snps> <path to gwas> <path to phenotype> <path to test> <test vec> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(pgenlibr)
}))

c_u = args[1] # causal uncorrected snps
a_u = args[2] # ascertained uncorrected snps
c_Tm = args[3] # causal Tm snps
a_Tm = args[4] # ascertained Tm snps
c_ID = args[5] # causal ID snps
a_ID = args[6] # ascertained ID snps
path_to_gwas = args[7] # path to GWAS genotypes
path_to_phenotype = args[8]
path_to_test = args[9]
tvec_file = args[10]
outfile = args[11]

# Read in SNPs
c_u <- fread(c_u)
a_u <- fread(a_u)
c_Tm <- fread(c_Tm)
a_Tm <- fread(a_Tm)
c_ID <- fread(c_ID)
c_ID <- fread(a_ID)

# Read in phenotypes
phenos <- fread(path_to_phenotype)

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

# Function to comput joint effect sizes
compute_joint <- function(geno_prefix, betas, phenos) {

  # Read in geotype matrix
  G <- read_genos(geno_prefix, betas)

  # Mean center genotype matrix
  G <- scale(G, scale = F)

  # Compute joint effect sizes
  mod <- lm(phenos$pheno_strat ~ G)
  betas_joint <- coef(mod)[-1]

  # Return joint effect sizes in a new column
  betas$joint <- betas_joint
  return(betas)
}

# Compute joint effect sizes
df_cu <- compute_joint(path_to_gwas, c_u, phenos)
df_au <- compute_joint(path_to_gwas, a_u, phenos)
df_cTm <- compute_joint(path_to_gwas, c_Tm, phenos)
df_aTm <- compute_joint(path_to_gwas, a_Tm, phenos)
df_cID <- compute_joint(path_to_gwas, c_ID, phenos)
df_aID <- compute_joint(path_to_gwas, a_ID, phenos)

# Load Test vector
std.tvec <- fread(tvec_file)
tvec <- std.tvec$Tvec

# Function to compute q
compute_q <- function(path_to_test, betas, tvec) {

  # Load Genotypes
  X <- read_genos(path_to_test, betas)

  # Mean center genotypes
  X <- scale(X, scale = F)

  # Compute q
  n <- nrow(X)
  L <- ncol(X)
  q_marginal <- (1/(n-1)) * (tvec %*% X %*% betas$BETA_Strat)
  q_joint <- (1/(n-1)) * (tvec %*% X %*% betas$joint)

  # Return
  out <- c(q_marginal, q_joint, L)
  return(out)
}


# Main function
main <- function(path_to_test, betas, tvec, type) {

  # Re-calculate q with different numbers of SNPs
  df <- as.data.frame(matrix(NA, ncol = 3, nrow = nrow(betas)))
  for (i in 1:nrow(betas)) {
    df[i, ] <- compute_q(path_to_test, betas[1:i, ], tvec)
  }
  df$type <- type

  return(df)
}

# Compute final output
out <- as.data.frame(matrix, ncol = 4, nrow= 6)
colnames(out) <- c("q_marginal", "q_joint", "L", "type")
out[1,] <- main(path_to_test, c_u, tvec, "uncorrected_causal")
out[2,] <- main(path_to_test, a_u, tvec, "uncorrected_ascertained")
out[3,] <- main(path_to_test, c_Tm, tvec, "Tm_causal")
out[4,] <- main(path_to_test, a_Tm, tvec, "Tm_ascertained")
out[5,] <- main(path_to_test, c_ID, tvec, "ID_causal")
out[6,] <- main(path_to_test, a_ID, tvec, "ID_ascertained")
fwrite(out, outfile,row.names=F,quote=F,sep="\t", col.names = T)





