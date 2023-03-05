library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("C", i)}
print(reps)
cases <- c("C1")
phenos <- c("LAT", "DIAG", "PS")
envs <- c("env-0.0", "env-0.2", "env-1.0")
dat <- expand.grid(reps, cases, phenos, envs)
colnames(dat) <- c("rep", "case", "pheno","env")

agg_all_data <- function(rep, dir_path, pheno, env, case) {
  
  gw <- fread(paste0(dir_path, rep,"/", case, "/h2-0/", pheno, "/",env, "/" , "pca_weights.txt"))
  gw$PC <- seq(1,nrow(gw))
  return(gw)

}

df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Simulate_Phenotypes/SimpleGrid/' )

fwrite(df, "SimpleGrid_pheno_weights.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)