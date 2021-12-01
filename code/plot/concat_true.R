library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("L", i)}
print(reps)
h2 <- c("same-0.3", "diff-0.3")
#envs <- c("env-0.0", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.05", "env-0.06", "env-0.07", "env-0.08", "env-0.09", "env-0.1")
#envs <- c("env-0.0","env-0.005", "env-0.01","env-0.015", "env-0.02","env-0.025", "env-0.03","env-0.035", "env-0.04","env-0.045", "env-0.05", "env-0.055","env-0.06")
envs <- c("env-0.0", "env-1.0")
cases <- c("C1")
ts <- c("p-0.50","p-0.55","p-0.60", "p-0.65", "p-0.70")
dat <- expand.grid(reps, cases, h2, envs, ts)
colnames(dat) <- c("rep", "case", "h2", "env", "ts")

agg_all_data <- function(rep, dir_path, case, type, h2, env,ts) {

  if (h2  == "diff-0.3" & env == "env-0.0" ) {
    Qx <- as.data.frame(matrix(NA, nrow=1,ncol=3))
    colnames(Qx) <- c("Qx", "P-Chi", "P-EN")
} else {  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", ts, "/",  env, "/", "Qx_true.txt"))
  colnames(Qx) <- c("Qx", "P-Chi", "P-EN")
}  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )

fwrite(df, "4PopSplit_L100_true.txt", row.names=F,quote=F,sep="\t", col.names = T)

#nc <- df %>% filter(type == "nc"| type == "nc-Tm" | type == "c" | type == "c-Tm")  %>% group_by(env, type, case) %>% summarise(fp_strat = sum(`P-EN` < 0.05)/ n(), #avg_Ax = mean(Ax))

#print(nc,n = Inf)