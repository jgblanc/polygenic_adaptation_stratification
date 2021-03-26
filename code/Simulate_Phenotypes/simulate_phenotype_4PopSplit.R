
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
env=as.numeric(args[5])
set.seed(args[6])

print(as.numeric(env))
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
prs$grandom = rnorm(sample_size,0, sqrt(1 - h2))

# effect on the "first" population (A)
pops <- unique(prs$pop)
prs$env = sapply(prs$pop,
                 function(x){
                   if(x == pops[1]){
                     rnorm(n = 1,
                           mean = env,
                           sd = sqrt(1 - h2)) }else{
                             rnorm(n = 1,
                                   mean = 0,
                                   sd = sqrt(1 - h2))
                           }})

#scale each so that their variances are 0.2 - to adjust heritability to h2
prs$grandom = scale( prs$grandom, scale = T) * sqrt( 1 - h2)
prs$env = scale(prs$env, scale = T) * sqrt(1 - h2)

#add prs to each of the environmental effects
prs = prs %>%
  mutate(pheno_random = prs + grandom,
        pheno_strat = prs + env)

fwrite(
  prs%>%
  mutate(FID=IID)%>%
  select(FID,IID,pheno_random,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
