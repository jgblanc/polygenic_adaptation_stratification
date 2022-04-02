library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "h2-0.0"
#envs <- c("env_0.0","env_0.005" ,"env_0.01","env_0.015","env_0.02","env_0.025", "env_0.03", "env_0.035","env_0.04", "env_0.045", "env_0.05", "env_0.055", "env_0.06")
envs <- c("env_0.0", "env_0.06")
cases <- c("C1", "C2")
ts <- c("p-0.50")
dat <- expand.grid(reps, cases, h2, ts,  envs)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/",  "Qx.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm", "c-ID", "c.p-ID", "nc-ID", "True")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "4PopSplit_h23_A.txt", row.names=F,quote=F,sep="\t", col.names = T)