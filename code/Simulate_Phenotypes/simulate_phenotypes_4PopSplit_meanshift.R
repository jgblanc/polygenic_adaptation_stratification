#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <hertaibility> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
popfile = args[2] # pop file from msprime simulation
output_file = args[3] #name of output file
h2 = as.numeric(args[4])
env_s=as.numeric(args[5])
set.seed(as.numeric(args[6]))

print(as.numeric(env_s))
print(as.numeric(h2))

prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

sample_size=nrow(prs)


#load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","pop")

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

##### simulate environmental effects
#no 'environmental' effect
#NOTE: SIGMA_G SET AT h2

# effect on the "first" population (A)
pops <- unique(prs$pop)
prs$env = rnorm(sample_size,0, sqrt(1 - h2))
prs <- prs %>% group_by(pop) %>% mutate(env = env-mean(env))
prs <- prs %>% group_by(pop) %>% mutate(env = ifelse(pop == pops[1], env + env_s, env ))

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
