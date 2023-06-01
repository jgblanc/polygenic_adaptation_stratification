library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "pc-0"
#envs <- c("env-0.0", "env-1.0", "env-0.2", "env-0.4", "env-0.6", "env-0.8",  "env-1.2", "env-1.4", "env-1.6", "env-1.8", "env-2.0")
envs  <-  c("env-2.0")
cases <- c("C1")
#test <- c("PS", "LAT")
test <- c("PS")
#pheno <- c("PS", "DIAG", "LAT")
pheno <- c("PS")
dat <- expand.grid(reps, cases, h2, pheno,  envs, test)
colnames(dat) <- c("rep", "case", "h2", "pheno", "envs", "test")

agg_all_data <- function(rep, dir_path, case, h2, pheno, envs, test) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", pheno, "/", envs, "/", test, "/",  "Qx.txt"))
  #Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm", "c-ID", "c.p-ID", "nc-ID")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/SimpleGrid/' )

fwrite(df, "SimpleGrid_fp_A_PCs.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)
