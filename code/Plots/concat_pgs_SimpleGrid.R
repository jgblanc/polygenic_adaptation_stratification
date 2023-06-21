library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "h2-0"
envs <- c("env-1.0", "env-0.2")
cases <- c("C1")
test <- c("LAT", "PS")
pheno <- c("LAT", "DIAG", "PS")
dat <- expand.grid(reps, cases, h2, pheno,  envs, test)
colnames(dat) <- c("rep", "case", "h2", "pheno", "envs", "test")

agg_all_data <- function(rep, dir_path, case, h2, pheno, envs, test) {
  
  PGS <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", pheno, "/", envs, "/", test, "/",  "PGS.txt"))
  pop <- fread(paste0('../../output/Simulate_Genotypes/SimpleGrid/', rep,"/", "genos.pop"))
  colnames(pop) <- c("IID", "#FID", "POP", "LAT", "LONG" )
  
  df <- inner_join(PGS, pop) %>% select("nc.joint", "nc_Tm.joint", "nc_10.joint", POP, LAT, LONG)
  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/SimpleGrid/' )

fwrite(df, "SimpleGrid_pgs_PCs_A.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)