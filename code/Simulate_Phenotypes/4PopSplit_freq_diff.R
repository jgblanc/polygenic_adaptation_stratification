args=commandArgs(TRUE)


suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
  library(tidyr)
}))

#frequency file
effects_file=args[1]
print(paste("The effects file is",effects_file))

#pop file
pop_file=args[2]
print(paste("The pop file is",pop_file))

#output file 
out_file = args[3]
print(paste("The output file is",out_file))

#test genotypes prefix
geno_prefix_test = args[4]

#gwas genotypes prefix
geno_prefix_gwas = args[5]

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

# Read in Betas
beta_df <- fread(effects_file)
colnames(beta_df) <- c("ID", "A1", "BETA")

# Read in population ID info
pop <- fread(pop_file, header = F)

# Get number of individuals in each test population
n1 <- as.numeric(count(pop,V3)[3,2])
n2 <- as.numeric(count(pop,V3)[4,2])

# Read in genotype matrix for causal variants
G <- read_genos(geno_prefix_test, beta_df[,"ID"])

# Calculate population specific allele frequeny
p1 <- colMeans(G[1:n1,])
p2 <- colMeans(G[(n1+1):(n1+n2),])

# Get allele frequency difference
diff_test <- p1 - p2

# Get number of individuals in each gwas population
n1 <- as.numeric(count(pop,V3)[1,2])
n2 <- as.numeric(count(pop,V3)[2,2])

# Read in genotype matrix for causal variants
G <- read_genos(geno_prefix_gwas, beta_df[,"ID"])

# Calculate population specific allele frequeny
p1 <- colMeans(G[1:n1,])
p2 <- colMeans(G[(n1+1):(n1+n2),])

# Get allele frequency difference
diff_gwas <- p1 - p2

# Create output table
df <- bind_cols(beta_df, diff_test, diff_gwas)
colnames(df) <- c("ID", "A1", "BETA", "Diff_Test", "Diff_GWAS")
fwrite(df, out_file,row.names=F,col.names=T,quote=F,sep="\t")






