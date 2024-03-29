library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "h2-0.3"
#envs <- c("env_0.0","env_0.0025","env_0.005", "env_0.0075","env_0.01","env_0.0125", "env_0.015","env_0.0175","env_0.02","env_0.0225", "env_0.025","env_0.0275","env_0.03", "env_0.0325", "env_0.035")
#envs <- c("env_0.0", "env_-0.1", "env_0.1", "env_0.2", "env_-0.2", "env_0.3", "env_-0.3")
envs <- c("env_0.0", "env_-0.1", "env_0.1")
cases <- c("C1")
ts <- c("p-0.50", "p-0.53", "p-0.56", "p-0.59", "p-0.62")
#ts <- c("p-0.50")
dat <- expand.grid(reps, cases, h2, ts,  envs)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs) {
  
  Q1 <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/",  "genos-gwas_common.mean"))
  Q1$type <- "uc"
  Q2 <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/",  "genos-gwas_common-Tm.mean"))
  Q2$type <- "Tm"
  Q3 <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/",  "genos-gwas_common-ID.mean"))
  Q3$type <- "ID"  

  Qx <- rbind(Q1, Q2, Q3)
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Run_GWAS/4PopSplit/' )

fwrite(df, "A_4PopSplit_ts_Bbins.txt", row.names=F,quote=F,sep="\t", col.names = T)
