#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <hertaibility> <environmental shift> <seed> <type of stratification>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
popfile = args[2] # pop file from msprime simulation
output_file = args[3] #name of output file
h2 = as.numeric(args[4]) # hertiability
env_s=as.numeric(args[5]) # magnitude of environmental shift
set.seed(as.numeric(args[6])) # random seed
pheno_pattern = args[7]


# Read in genetic values
prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

# Get simulation sample size
sample_size=nrow(prs)

# Set larger GWAS sample size (we will match errors to this sample size)
N_GWAS = 500000 #### Fix to be per pop

# Load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","Pop", "Lat", "Long")

# Add this info to genetic value file
prs=merge(prs, pop, by="IID", sort=F)

# Calculate average phenotype per pop
prs$env = rnorm(sample_size,0, sqrt(1 - h2))
Z <- prs %>% group_by(Pop) %>% summarise(avg = mean(env)) %>% pull(avg)

# Rescale averages to larger GWAS size
n_sim <- prs %>% group_by(Pop) %>% summarise(num = n()) %>% pull(num)
Z_gwas <- sqrt(n_sim/ ((n_sim/sum(n_sim) * N_GWAS))) * Z

# Rescale individual phenotypes
delta <- Z - Z_gwas
pops <- seq(0, length(delta)-1)
for (i in 1:length(delta)) {
  prs = prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == pops[i], env - delta[i], env))
}

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

# Add genetic value to each of the environmental effects
prs = prs %>%
  mutate(pheno_strat = prs + env)

fwrite(
  prs%>%
    mutate(FID=IID)%>%
    ungroup() %>%
    select(FID,IID,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
