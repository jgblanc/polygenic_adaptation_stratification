# This script takes in selected SNPs and re-estimates effect sizes jointly

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript joint_effect_sizes.R <path to gwas genotype> <path to SNPs>
                        <paht to phenotype> <path to covariates> <list of PCs>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(pgenlibr)
}))


path_to_gwas = args[1] # path to GWAS genotypes
path_to_SNPs = args[2]
path_to_phenotype = args[3]
covar_file = args[4]
pc_list = args[5]
print(pc_list)
pcs <- as.numeric(strsplit(pc_list,"-")[[1]])

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
compute_joint <- function(geno_prefix, betas, phenos, covar) {

  # Read in geotype matrix
  G <- read_genos(geno_prefix, betas)

  # Mean center genotype matrix
  G <- scale(G, scale = F)

  # Compute joint effect sizes
  mod <- lm(phenos$pheno_strat ~ G + covar)
  betas_joint <- coef(mod)[2:(ncol(G)+1)]

  # Compute marginal
  betas_marginal <- rep(0, nrow(betas))
  for (i in 1:nrow(betas)){
    betas_marginal[i] <- coef(lm(phenos$pheno_strat ~ G[,i] + covar))[2]
  }

  # Return joint effect sizes in a new column
  betas$joint <- betas_joint
  betas$marginal <- betas_marginal
  return(betas)
}

# Compute joint effect sizes and save
df_covar <- fread(covar_file)

## Uncorrected
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, ".c.betas")), phenos, rep(0, nrow(phenos)))
fwrite(df, paste0(path_to_SNPs, ".c.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, ".nc.betas")), phenos, rep(0, nrow(phenos)))
fwrite(df, paste0(path_to_SNPs, ".nc.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
print("Finished uncorrected")

## Tm corrected
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-Tm.c.betas")), phenos, df_covar$Tm)
fwrite(df, paste0(path_to_SNPs, "-Tm.c.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-Tm.nc.betas")), phenos, df_covar$Tm)
fwrite(df, paste0(path_to_SNPs, "-Tm.nc.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
print("Finished Tm")

## ID corrected
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-ID.c.betas")), phenos, df_covar$PopID)
fwrite(df, paste0(path_to_SNPs, "-ID.c.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-ID.nc.betas")), phenos, df_covar$PopID)
fwrite(df, paste0(path_to_SNPs, "-ID.nc.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
print("Finished ID")

## PC corrected
for (i in pcs) {

  df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-",i,".c.betas")), phenos, as.matrix(df_covar[,c(3,5:(i+4))]))
  fwrite(df, paste0(path_to_SNPs, "-", i,"_pcFGr", ".c.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)
  df <- compute_joint(path_to_gwas, fread(paste0(path_to_SNPs, "-", i, ".nc.betas")), phenos, as.matrix(df_covar[,c(3,5:(i+4))]))
  fwrite(df, paste0(path_to_SNPs, "-", i,"_pcFGr", ".nc.betas.joint"),row.names=F,quote=F,sep="\t", col.names = T)

}
print("Finished PC")







