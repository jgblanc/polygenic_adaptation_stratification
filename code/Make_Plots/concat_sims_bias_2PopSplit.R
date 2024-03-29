library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "h2-0.0"
#envs <- c("env_0.0", "env_0.01", "env_0.02", "env_0.03", "env_0.04", "env_0.05", "env_0.06", "env_0.07", "env_0.08","env_0.09", "env_0.1")
envs <- c("env_0.0","env_0.01", "env_0.02", "env_0.03", "env_0.04", "env_0.05", "env_0.06", "env_0.07", "env_0.08", "env_0.09", "env_0.1")
cases <- c("C1", "C2")
ts <- c("p-0.50")
nc <- c("c-200","c-100")
dat <- expand.grid(reps, cases, h2, ts,  envs, nc)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs", "nc")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs, nc) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", nc, "/", envs, "/",  "Qx.txt"))
  #Qx$type <- c("nc-uncorrected", "nc-Tm", "nc-ID", "c-uncorrected", "c-Tm","c-ID","true")

  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "A_2PopSplit_bias.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)
