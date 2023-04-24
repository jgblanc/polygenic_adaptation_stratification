# Split individuals into GWAS and test panel

args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


ss_test = as.numeric(args[1])
input_file = args[2]
gwas_file = args[3]
test_file = args[4]

# Read in pop info
pop <- fread(input_file, header = F)
colnames(pop) <- c("#IID", "FID", "POP", "LAT", "LONG" )

# Sample ss_test individuals per deme
df_test <- pop %>% group_by(POP) %>% sample_n(ss_test)
test <- df_test[,1:2]

# Use the rest of the individuals for the GWAS panel
df_gwas <- anti_join(pop, df_test)
gwas <- df_gwas[,1:2]

# Write panel IDs to file
write.table(gwas, gwas_file, quote = F, col.names = F, row.names = F)
write.table(test, test_file, quote = F, col.names = F, row.names = F)

