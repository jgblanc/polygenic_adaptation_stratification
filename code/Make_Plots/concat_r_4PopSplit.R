library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
cases <- c("C1", "C2")
dat <- expand.grid(reps, cases)
colnames(dat) <- c("rep", "case")

agg_all_data <- function(rep, dir_path, case) {
  
  Qx <- fread(paste0(dir_path, rep,"/", case, "/",  "r_bins.txt"))  

  return(Qx)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Calculate_Tm/4PopSplit/' )

fwrite(df, "A_4PopSplit_rbins.txt", row.names=F,quote=F,sep="\t", col.names = T)
