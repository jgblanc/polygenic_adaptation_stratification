library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("C", i)}
print(reps)
h2 <- "shift-0"
envs <- c("env-0.0", "env-0.5", "env-1.0")
cases <- c("C1")
pheno <- c("LAT", "DIAG", "PS")
dat <- expand.grid(reps, cases, h2, pheno,  envs)
colnames(dat) <- c("rep", "case", "h2", "pheno", "envs")

agg_all_data <- function(rep, dir_path, case, h2, pheno, envs) {
  
  phenos <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", pheno, "/", envs, "/", "genos-gwas_common.phenos.txt"))
  colnames(phenos) <- c("#FID", "IID", "pheno_strat")
  pop <- fread(paste0('../../output/Simulate_Genotypes/SimpleGrid/', rep,"/", "genos.pop"))
  colnames(pop) <- c("IID", "#FID", "POP", "LAT", "LONG" )
  
  df <- inner_join(phenos, pop)%>% select("pheno_strat", POP, LAT, LONG)
  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Simulate_Phenotypes/SimpleGrid/' )

fwrite(df, "SimpleGrid_phenos_C_shift.txt", row.names=F,quote=F,sep="\t", col.names = T)