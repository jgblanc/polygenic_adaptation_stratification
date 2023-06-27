library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
cases <- c("C1")
test <- c("LAT", "PS")
dat <- expand.grid(reps, cases, test)
colnames(dat) <- c("rep", "case", "test")

agg_all_data <- function(rep, dir_path, case, test) {
  
  gw <- fread(paste0(dir_path, rep,"/", case, "/", test, "/",  "gwas_pca_weights.txt"))
  gw$PC <- seq(1,nrow(gw))
  return(gw)

}

df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Calculate_FGr/SimpleGrid/' )

fwrite(df, "SimpleGrid_gw_scaleEvecs_A.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)