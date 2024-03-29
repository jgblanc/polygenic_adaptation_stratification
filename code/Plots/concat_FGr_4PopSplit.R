library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("B", i)}
print(reps)
cases <- c("C1", "C2") 
dat <- expand.grid(reps, cases)
colnames(dat) <- c("rep", "case")

agg_all_data <- function(rep, dir_path, case) {
  
  FGr <- fread(paste0(dir_path, rep,"/", case, "/", "Tm.txt"))
  IDs <- fread(paste0("../../output/Calculate_FGr/4PopSplit/B1/", "/", case, "/", "Tm-ID_covars.txt"))
  IDs <- IDs[,1:2]
  IDs$Tm <- FGr$Tm

  pop <- fread(paste0('../../output/Simulate_Genotypes/4PopSplit/', rep,"/", "genos.pop"))
  colnames(pop) <- c("IID", "#FID", "PopID" )
  
  df <- inner_join(pop, IDs)    

  return(df)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Calculate_FGr/4PopSplit/' )

fwrite(df, "4PopSplit_FGr_B.txt", row.names=F,quote=F,sep="\t", col.names = T)