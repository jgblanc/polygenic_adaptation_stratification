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
a_ID <- fread(a_ID)

# Read in phenotypes
phenos <- fread(path_to_phenotype)

print("Loaded data")

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
#df_cu <- compute_joint(path_to_gwas, c_u, phenos)
#df_au <- compute_joint(path_to_gwas, a_u, phenos)
#df_cTm <- compute_joint(path_to_gwas, c_Tm, phenos)
#df_aTm <- compute_joint(path_to_gwas, a_Tm, phenos)
#df_cID <- compute_joint(path_to_gwas, c_ID, phenos)
#df_aID <- compute_joint(path_to_gwas, a_ID, phenos)
#print("Computed joint effects")

# Load Test vector
std.tvec <- fread(tvec_file)
tvec <- std.tvec$Tvec

# Function to compute q
compute_q <- function(path_to_test, betas, tvec) {

  # Load Genotypes
  X <- read_genos(path_to_test, betas)
  freq <- colMeans(X)/2

  # Mean center genotypes
  X <- scale(X, scale = F)

  # Compute q
  n <- nrow(X)
  L <- ncol(X)

  q_marginal <- (1/(n-1)) * (tvec %*% X %*% as.matrix(betas$BETA_Strat))
  q_joint <- (1/(n-1)) * (tvec %*% X %*% as.matrix(betas$joint))

  # Compute Va
  Va_marginal <- 4 * sum(betas$BETA_Strat^2 * freq  * (1 - freq))
  Va_joint <- 4 * sum(betas$joint^2 * freq  * (1 - freq))

  # Return
  out <- c(q_marginal, q_joint,  Va_marginal, Va_joint, L)
  return(out)
}


# Main function
main <- function(path_to_test, betas, tvec, type) {

  # Re-calculate q with different numbers of SNPs
  df <- as.data.frame(matrix(NA, ncol = 5, nrow = nrow(betas)))
  for (i in 1:nrow(betas)) {

    # Compute joint effects
    b <- compute_joint(path_to_gwas, betas[1:i, ], phenos)

    # Compute q
    df[i, ] <- compute_q(path_to_test, b, tvec)
  }
  df$type <- type

  return(df)
}



# Compute final output
out <- main(path_to_test, a_u , tvec, "uncorrected_ascertained")
#out <- rbind(out, main(path_to_test, c_u, tvec, "uncorrected_ascertained"))
#out <- rbind(out, main(path_to_test, df_cTm, tvec, "Tm_causal"))
#out <- rbind(out, main(path_to_test, df_aTm, tvec, "Tm_ascertained"))
#out <- rbind(out, main(path_to_test, df_cID, tvec, "ID_causal"))
#out <- rbind(out, main(path_to_test, df_aID, tvec, "ID_ascertained"))
colnames(out) <- c("q_marginal", "q_joint", "Va_marginal", "Va_joint", "L", "type")
print(out)
fwrite(out, outfile,row.names=F,quote=F,sep="\t", col.names = T)





