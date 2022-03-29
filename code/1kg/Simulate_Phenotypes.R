#script takes in genetic values and generates phenotypes with and without stratification

args=commandArgs(TRUE)

if(length(args)<6){stop("Rscript simphenotype_ge.R <genetic value file> <pop file> <output_file> <hertaibility> <environmental shift> <eigenvectors> <eigenvalues>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
popfile = args[2] # pop file from msprime simulation
output_file = args[3] #name of output file
h2 = as.numeric(args[4]) #target heritability
evec_file = args[5]
eval_file = args[6]

print(paste0("The heritability is " , as.numeric(h2)))

# Load true genetic values
prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

#load file containing the pop ID for each inidividual
pop=fread(popfile,header=T)
colnames(pop)<-c("IID","pop")

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

# Draw random environment
sample_size <- nrow(prs)
prs$env = rnorm(sample_size,0, sqrt(1 - h2))
prs$env = scale(prs$env, scale = TRUE) * (sqrt(1 - h2))

# Load eigenvectors and values
vecs <- fread(evec_file)
vec_mat <- apply(vecs[,2:ncol(vecs)], 2, as.numeric)
vals <- fread(eval_file)
weighted_vecs <- vec_mat %*% diag(vals$V1)

# Use weighted sum of eigenvectors to add stratification
#pc_w <- rbeta(nrow(prs) -1, 1, 3)
pc_w <- c(1, rep(0,nrow(prs) -2 ))
strat_env <- weighted_vecs %*% pc_w
prs$strat <- strat_env[,1]
prs$pheno_strat <- prs$prs + prs$env + prs$strat


# Print heritability
print(paste("h2 :",
            round(cor(prs$prs,
                      prs$pheno_strat)^2,2)))

# Print mean phenotype in each pop
#print(prs %>% group_by(pop) %>% summarise(mean(pheno_strat)))

# Write to output file
fwrite(
  prs%>%
    mutate(FID=IID)%>%
    ungroup() %>%
    select(FID,IID,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
