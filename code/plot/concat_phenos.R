library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("T", i)}
print(reps)
h2 <- "cell-0"
#envs <- c("env-0.0", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.05", "env-0.06", "env-0.07", "env-0.08", "env-0.09", "env-0.1")
envs <- "env-0.0"
cases <- c("C1")
dat <- expand.grid(reps, cases, h2, envs)
colnames(dat) <- c("rep", "case", "h2", "env")

agg_all_data <- function(rep, dir_path, case, type, h2, env) {
  
  pheno <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", env, "/", "genos-gwas_common.phenos.txt"))
  pop <- fread(paste0('../../output/Simulate_Genotypes/SimpleGrid/', rep,"/", "genos.pop"))
  colnames(pop) <- c("IID", "FID", "POP", "LAT", "LONG" )
  
  df <- inner_join(pheno, pop)%>% select("IID", "pheno_strat", "POP", "LAT", "LONG")
  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Simulate_Phenotypes/SimpleGrid/' )


fwrite(df, "SimpleGrid_cell_T100_phenos.txt", row.names=F,quote=F,sep="\t", col.names = T)