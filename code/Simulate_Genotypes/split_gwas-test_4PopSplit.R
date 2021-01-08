# Split Population IDs into two sets of test/gwas pairs

print("Running")
library(data.table)

# Read in pop info
pop <- fread(snakemake@input[[1]], header = F)

# Case 1: A,C gwas B,D test
gwas <- subset(pop, pop$V3 == "A" | pop$V3 == "C")
test <- subset(pop, pop$V3 == "B" | pop$V3 == "D")

write.table(gwas, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
write.table(test, snakemake@output[[2]], quote = F, col.names = F, row.names = F)

# Case 2: A,B gwas, C,D test
gwas <- subset(pop, pop$V3 == "A" | pop$V3 == "B")
test <- subset(pop, pop$V3 == "C" | pop$V3 == "D")

write.table(gwas, snakemake@output[[3]], quote = F, col.names = F, row.names = F)
write.table(test, snakemake@output[[4]], quote = F, col.names = F, row.names = F)
