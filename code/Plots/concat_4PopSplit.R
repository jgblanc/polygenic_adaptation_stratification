library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("B", i)}
print(reps)
h2 <- "h2-0.3"
#envs <- c("env_0.0", "env_0.01", "env_0.02", "env_0.03", "env_0.04", "env_0.05", "env_0.06", "env_0.07", "env_0.08","env_0.09", "env_0.1")
envs <- c("env_0.0","env_-0.2", "env_0.2")
cases <- c("C1")
ts <- c("p-0.50", "p-0.53", "p-0.56", "p-0.59", "p-0.62")
#ts <- "p-0.50"
nc <- c("c-200")
dat <- expand.grid(reps, cases, h2, ts,  envs, nc)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs", "nc")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs, nc) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", nc, "/", envs, "/",  "Qx_ss.txt"))
  #Qx$type <- c("nc-uncorrected", "nc-Tm", "nc-ID", "c-uncorrected", "c-Tm","c-ID","true")

  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "ts_ss.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)
