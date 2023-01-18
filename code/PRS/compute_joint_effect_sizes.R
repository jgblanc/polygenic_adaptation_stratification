# This script takes in selected SNPs, re-estimates effect sizes jointly, and computes q for different numbers of SNPs

args=commandArgs(TRUE)

if(length(args)<9){stop("Rscript joint_effect_sizes.R <asecertained snps> <causal snps> <path to gwas> <path to phenotype> <path to test> <test vec> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(pgenlibr)
}))

c_u_file = args[1] # causal uncorrected snps
a_u_file = args[2] # ascertained uncorrected snps
c_Tm_file = args[3] # causal Tm snps
a_Tm_file = args[4] # ascertained Tm snps
c_ID_file = args[5] # causal ID snps
a_ID_file = args[6] # ascertained ID snps
path_to_gwas = args[7] # path to GWAS genotypes
path_to_phenotype = args[8]
path_to_test = args[9]


# Read in SNPs
c_u <- fread(c_u_file)
a_u <- fread(a_u_file)
c_Tm <- fread(c_Tm_file)
a_Tm <- fread(a_Tm_file)
c_ID <- fread(c_ID_file)
a_ID <- fread(a_ID_file)

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

  # Compute marginal
  betas_marginal <- rep(0, nrow(betas))
  for (i in 1:nrow(betas)){
    betas_marginal[i] <- coef(lm(phenos$pheno_strat ~ G[,i]))[2]
  }

  # Return joint effect sizes in a new column
  betas$joint <- betas_joint
  betas$marginal <- betas_marginal
  return(betas)
}

# Compute joint effect sizes
df_cu <- compute_joint(path_to_gwas, c_u, phenos)
df_au <- compute_joint(path_to_gwas, a_u, phenos)
df_cTm <- compute_joint(path_to_gwas, c_Tm, phenos)
df_aTm <- compute_joint(path_to_gwas, a_Tm, phenos)
df_cID <- compute_joint(path_to_gwas, c_ID, phenos)
df_aID <- compute_joint(path_to_gwas, a_ID, phenos)
print("Computed joint effects")

# Save output
fwrite(df_cu, paste0(c_u_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_au, paste0(a_u_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_cTm, paste0(c_Tm_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_aTm, paste0(a_Tm_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_cID, paste0(c_ID_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)
fwrite(df_aID, paste0(a_ID_file, ".joint"),row.names=F,quote=F,sep="\t", col.names = T)






