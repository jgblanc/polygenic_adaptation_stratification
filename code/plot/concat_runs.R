library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("F", i)}
#reps = reps[!(reps %in% c(5))]
print(reps)
h2 <- "h2-0"
#envs <- c("env-0.0", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.05", "env-0.06", "env-0.07", "env-0.08", "env-0.09", "env-0.10")
envs <- c("env-0.0", "env-0.5")
cases <- c("C1")
dat <- expand.grid(reps, cases, h2, envs)
colnames(dat) <- c("rep", "case", "h2", "env")

agg_all_data <- function(rep, dir_path, case, type, h2, env) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", env, "/", "Qx.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/SimpleGrid/' )

fwrite(df, "F100_SimpleGrid.txt", row.names=F,quote=F,sep="\t", col.names = T)