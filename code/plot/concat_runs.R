library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("T", i)}
#reps = reps[!(reps %in% c(5))]
print(reps)
h2 <- "3Pop-0"
#envs <- c("env-0.0", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.05", "env-0.06", "env-0.07", "env-0.08", "env-0.09", "env-0.1")
#envs <- c("env-0.0","env-0.005", "env-0.01","env-0.015", "env-0.02","env-0.025", "env-0.03","env-0.035", "env-0.04","env-0.045", "env-0.05", "env-0.055","env-0.06")
envs <- c('env-0.0', "env-0.3")
cases <- c("C2")
dat <- expand.grid(reps, cases, h2, envs)
colnames(dat) <- c("rep", "case", "h2", "env")

agg_all_data <- function(rep, dir_path, case, type, h2, env) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/", h2, "/", env, "/", "Qx_one.txt"))
  Qx$type <- c("c", "c.p", "nc", "c-Tm", "c.p-Tm", "nc-Tm")
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/SimpleGrid/' )

fwrite(df, "SimpleGrid_T100_C2_one.txt", row.names=F,quote=F,sep="\t", col.names = T)

nc <- df %>% filter(type == "nc"| type == "nc-Tm") %>% filter(case == "C2")  %>% group_by(env, type, case) %>% summarise(fp_strat = sum(`P-EN` < 0.05)/ n(), avg_Ax = mean(Ax))

print(nc)