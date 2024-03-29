args=commandArgs(TRUE)

# Split Population IDs into two sets of test/gwas pairs

library(data.table)
library(dplyr)

set.seed(as.numeric(args[7]))


# Read in pop info
pop <- fread(args[1], header = T)

# Get test set size
size = as.numeric(args[2])
print(size)

# Case 1: A,C
gwas <- subset(pop, pop$POP == "A" | pop$POP == "C")
test <- gwas %>% group_by(POP) %>% sample_n(size/2)

# Remove test inds from gwas
gwas <- gwas  %>% filter(!FID %in% test$FID)

write.table(gwas, args[3], quote = F, col.names = T, row.names = F)
write.table(test, args[4], quote = F, col.names = T, row.names = F)

# Case 2: A,B gwas, C,D test
gwas <- subset(pop, pop$POP == "A" | pop$POP == "B")
test <- gwas %>% group_by(POP) %>% sample_n(size/2)

# Remove test inds from gwas
gwas <- gwas  %>% filter(!FID %in% test$FID)

write.table(gwas, args[5], quote = F, col.names = T, row.names = F)
write.table(test, args[6], quote = F, col.names = T, row.names = F)
