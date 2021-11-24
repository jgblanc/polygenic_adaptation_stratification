library(dplyr)
library(data.table)
print("Running")

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("T", i)}
print(reps)
h2 <- "same_large-0.3"
envs <- c("env-0.0","env-1.0")
cases <- c("C1")
ts <- c("p-0.50","p-0.55","p-0.60", "p-0.65", "p-0.70")
dat <- expand.grid(reps, cases, h2, envs, ts)
colnames(dat) <- c("rep", "case", "h2", "env", "ts")
print(dat)

agg_all_data <- function(rep, dir_path, case, type, h2, env, ts) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/", env, "/", "Qx_ID.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm", "c-ID", "c.p-ID", "nc-ID")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "4PopSplit_T100_true_same_large.txt", row.names=F,quote=F,sep="\t", col.names = T)

#nc <- df %>% filter(type == "nc"| type == "nc-Tm" | type == "c" | type == "c-Tm")  %>% group_by(env, type, case, ts) %>% summarise(fp_strat = sum(`P-EN` < 0.05)/ n(), avg_Ax = mean(Ax))

#print(nc,n = Inf)