#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <hertaibility> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
#gvalue_file = "output/Simulate_Phenotypes/SimpleGrid/E2/C1/h2-0/genos-gwas_common.gvalue.sscore"
popfile = args[2] # pop file from msprime simulation
#popfile = "output/Simulate_Genotypes/SimpleGrid/E2/genos.pop"
output_file = args[3] #name of output file
h2 = as.numeric(args[4])
env_s=as.numeric(args[5])
print(env_s)
set.seed(as.numeric(args[6]))


prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

sample_size=nrow(prs)
N_GWAS = 500000

#load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","Pop", "Lat", "Long")

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

##### simulate environmental effects
#no 'environmental' effect
#NOTE: SIGMA_G SET AT h2

# Shift the environmental variable by env_s/lat
num_lat_demes <- length(unique(prs$Lat))
if (env_s == 0) {
  shifts <- rep(0, num_lat_demes)
} else {
  shifts <- seq(0, env_s, env_s / (num_lat_demes-1))
}
print(shifts)

# Shift environment for individuals in each deme
#prs$env = rnorm(sample_size,0, sqrt(1 - h2))
#prs <- prs %>% group_by(Lat) %>% mutate(env = env-mean(env))
#prs <- prs %>% group_by(Pop) %>% mutate(env = env-mean(env))
#for (i in 1:num_lat_demes) {
#  prs <- prs %>% group_by(Lat) %>% mutate(env = ifelse(Lat == (i-1), env + shifts[i], env))
#}

# Calculate average phenotype per pop
prs$env = rnorm(sample_size,0, sqrt(1 - h2))
Z <- prs %>% group_by(Pop) %>% summarise(avg = mean(env)) %>% pull(avg)

# Rescale to GWAS size
n_sim <- prs %>% group_by(Pop) %>% summarise(num = n()) %>% pull(num)
Z_gwas <- sqrt(n_sim/ ((n_sim/sum(n_sim) * N_GWAS))) * Z

# Rescale individual phenotypes
delta <- Z - Z_gwas
pops <- seq(0, length(delta)-1)
for (i in 1:length(delta)) {
  prs = prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == pops[i], env - delta[i], env))
}

# Add stratification along diagnoal Latitude
id_diag <- c(0,7,14, 21, 28, 35)
for (i in 1:length(id_lat)) {
  prs <- prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == id_diag[i], env + shifts[i], env))
}

#add prs to each of the environmental effects
prs = prs %>%
  mutate(pheno_strat = prs + env)

fwrite(
  prs%>%
    mutate(FID=IID)%>%
    ungroup() %>%
    select(FID,IID,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
