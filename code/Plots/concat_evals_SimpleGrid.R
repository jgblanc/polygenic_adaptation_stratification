library(dplyr)
library(data.table)

n <- 100
reps <- rep(NA, n)
for (i in 1:n){reps[i] <- paste0("A", i)}
print(reps)
cases <- c("C1") 
dat <- expand.grid(reps, cases)
colnames(dat) <- c("rep", "case")

agg_all_data <- function(rep, dir_path, case) {
  
  evals <- fread(paste0(dir_path, rep,"/", case, "/", "genos-gwas.eigenval"))
  evals$number  <- seq(1, nrow(evals))
  return(evals)
}
df <- plyr::mdply(dat, agg_all_data, dir_path = '../../output/Calculate_FGr/SimpleGrid/' )

fwrite(df, "SimpleGrid_evals_A.txt", row.names=F,quote=F,sep="\t", col.names = T)