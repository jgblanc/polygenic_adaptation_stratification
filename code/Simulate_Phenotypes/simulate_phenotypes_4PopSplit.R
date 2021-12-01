#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <hertaibility> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
#gvalue_file = 'output/Simulate_Phenotypes/4PopSplit/E2/C1/h2-0/genos-gwas_common.gvalue.sscore'
popfile = args[2] # pop file from msprime simulation
#popfile = "output/Simulate_Genotypes/4PopSplit/E2/genos.pop"
output_file = args[3] #name of output file
h2 = as.numeric(args[4])
env_s=as.numeric(args[5])
#set.seed(as.numeric(args[6]))

print(as.numeric(env_s))
print(as.numeric(h2))

prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

sample_size=nrow(prs)
NGWAS = 500000
N1_gwas = NGWAS / 2
N2_gwas = NGWAS / 2

#load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","pop")

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

##### simulate environmental effects
#no 'environmental' effect
#NOTE: SIGMA_G SET AT h2

# effect on the "first" population (A)
#pops <- unique(prs$pop)
#prs$env = rnorm(sample_size,0, sqrt(1 - h2))
#prs <- prs %>% group_by(pop) %>% mutate(env = env-mean(env))
#prs <- prs %>% group_by(pop) %>% mutate(env = ifelse(pop == pops[1], env + env_s, env ))

# Draw random environment
pops <- unique(prs$pop)
prs$env = rnorm(sample_size,0, sqrt(1 - h2))
prs$env = scale(prs$env, scale = TRUE) * (sqrt(1 - h2))

# Add stratification effect to environment
prs <- prs %>% group_by(pop) %>% mutate(env = ifelse(pop == pops[1], env + env_s, env ))

# Calculate average phenotype
Z1 = prs %>% group_by(pop) %>% summarise(avg = mean(env)) %>% filter(pop == pops[1]) %>% pull(avg)
Z2 = prs %>% group_by(pop) %>% summarise(avg = mean(env)) %>% filter(pop == pops[2]) %>% pull(avg)

# Rescale averages
N1 = nrow(prs %>% filter(pop == pops[1]))
N2 = nrow(prs %>% filter(pop == pops[2]))
Z1_gwas = sqrt(N1/N1_gwas) * Z1
Z2_gwas = sqrt(N2/N2_gwas) * Z2

# Recalulate individual phenotypes
delta1 = Z1 - Z1_gwas
delta2 = Z2 - Z2_gwas
prs = prs %>% group_by(pop) %>% mutate(env = ifelse(pop == pops[1], env - delta1, env - delta2))

# Add environmental effect to
prs = prs %>%
  mutate(pheno_strat = prs + env)


print(paste("h2 :",
            round(cor(prs$prs,
                      prs$pheno_strat)^2,2)))
fwrite(
  prs%>%
    mutate(FID=IID)%>%
    ungroup() %>%
    select(FID,IID,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
