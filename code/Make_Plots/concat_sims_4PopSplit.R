library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
h2 <- "h2-0.3"
#envs <- c("env_0.0","env_0.0025","env_0.005", "env_0.0075","env_0.01","env_0.0125", "env_0.015","env_0.0175","env_0.02","env_0.0225", "env_0.025","env_0.0275","env_0.03", "env_0.0325", "env_0.035")
#envs <- c("env_0.0","env_0.01","env_0.02","env_0.03","env_0.04", "env_0.05","env_0.06","env_0.07","env_0.08","env_0.09", "env_0.1", "env_0.11","env_0.12", "env_0.13", "env_0.14", "env_0.15")
envs <- c("env_0.0", "env_-0.1", "env_0.1", "env_0.2", "env_-0.2", "env_0.3", "env_-0.3")
#envs <- c("env_0.0", "env_1.0", "env_-1.0")
cases <- c("C1", "C2")
#ts <- c("p-0.50", "p-0.53", "p-0.56", "p-0.59", "p-0.62")
ts <- c("p-0.50")
dat <- expand.grid(reps, cases, h2, ts,  envs)
colnames(dat) <- c("rep", "case", "h2", "ts", "envs")

agg_all_data <- function(rep, dir_path, case, h2, ts, envs) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", envs, "/",  "Qx.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm", "c-ID", "c.p-ID", "nc-ID", "True")
  #Qx$type <- c("Uncorrected", "Tm", "ID")  

  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "A_4PopSplit_h23.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)
