library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("D", i)}
print(reps)
dat <- expand.grid(reps)
colnames(dat) <- c("rep")

agg_all_data <- function(rep, dir_path) {
  
  Qx <- fread(paste0(dir_path, rep,"/C1/", "q.txt"))
  #print(head(Qx))
  
  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/PGA_test/4PopSplit/' )
#colnames(df) <- c("rep", "case", "num","true.q", "S.q", "L.q","ncausal", "eq")

fwrite(df, "Sig_Results_Indep.txt.gz", row.names=F,quote=F,sep="\t", col.names = T)
