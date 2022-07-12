library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "snps-0.0"
envs <- c("env_0.03")
cases <- c("C1")
ts <- c("p-0.50")
snps <-c("snp-500")
dat <- expand.grid(reps, cases, h2, ts,  envs, snps)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs", "snps")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs, snps) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/", snps, "/", "Qx.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm", "c-ID", "c.p-ID", "nc-ID", "True")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "A_4PopSplit_snps.txt", row.names=F,quote=F,sep="\t", col.names = T)
