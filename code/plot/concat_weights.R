library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("T", i)}
print(reps)
h2 <- "scale-0"
#envs <- c("env-0.0", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.05", "env-0.06", "env-0.07", "env-0.08", "env-0.09", "env-0.1") 
cases <- c("C1")
dat <- expand.grid(reps, cases)
colnames(dat) <- c("rep", "case")

agg_all_data <- function(rep, dir_path, case) {
  
  w <- fread(paste0(dir_path, rep,"/", case, "/", "gwas_pca_weights.txt"))
  w$PC <- seq(1, nrow(w))  
  return(w)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Calculate_Tm/SimpleGrid/' )


fwrite(df, "SimpleGrid_LatLat_weights.txt", row.names=F,quote=F,sep="\t", col.names = T)