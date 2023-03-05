## This script calculates the relative weights of GWAS PCs on T^GWAS

## Requires: gwas_pca.eigenvec, Tm.txt

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_Tm.R <eigenvecs> <Tm> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))


vecs_file = args[1] # eigenvectors
pheno_type = args[2] # Tm
popfile = args[3]
gvalue_file = args[4]
out_file = args[5]

# Load gwas eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)

## Get confounder

# Read in genetic values
prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

# Get simulation sample size
sample_size=nrow(prs)


# Load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","Pop", "Lat", "Long")

# Add this info to genetic value file
prs=merge(prs, pop, by="IID", sort=F)
prs$env = 0

if (pheno_pattern == "LAT") {

  print(pheno_pattern)

  # Shift the environmental variable by env_s/lat
  num_lat_demes <- length(unique(prs$Lat))
  if (env_s == 0) {
    shifts <- rep(0, num_lat_demes)
  } else {
    shifts <- seq(0, env_s, env_s / (num_lat_demes-1))
  }

  # Add stratification along Latitude
  for (i in 1:num_lat_demes) {
    prs <- prs %>% group_by(Lat) %>% mutate(env = ifelse(Lat == (i-1), env + shifts[i], env))
  }

} else if (pheno_pattern == "DIAG") {

  print(pheno_pattern)

  # Shift the environmental variable by env_s/# of diagonal squares
  num_diag_demes <- length(unique(prs$Lat))
  if (env_s == 0) {
    shifts <- rep(0, num_diag_demes)
  } else {
    shifts <- seq(0, env_s, env_s / (num_diag_demes-1))
  }

  # Add stratification along diagonal Latitude
  id_diag <- c(0,7,14, 21, 28, 35)
  for (i in 1:length(id_diag)) {
    prs <- prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == id_diag[i], env + shifts[i], env))
  }
} else if (pheno_pattern == "PS") {

  print(pheno_pattern)

  # Add stratification on deme 25
  id_diag <- c(25)
  for (i in 1:length(id_diag)) {
    prs <- prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == id_diag[i], env + env_s, env))
  }

} else {
  stop("Please enter acceptable phenotype pattern: LAT, DIAG, PS")
}

Tm = prs$env

# Get weights
B = t(vecs) %*% Tm

# Get relative weights
Br = abs(B) / sum(abs(B))

# Convert to dataframe
B <- as_tibble(cbind(B, Br))
colnames(B) <- c("PC_weights", "Relative_PC_weights")

# Save weights
fwrite(B, out_file,row.names=F,quote=F,sep="\t", col.names = T)
