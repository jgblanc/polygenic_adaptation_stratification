# Split individuals into GWAS and test panel

args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


ss_test = as.numeric(args[1])
input_file = args[2]
gwas_north_file = args[3]
gwas_south_file = args[4]
gwas_joint_file = args[5]
test_file = args[6]

# Read in pop info
pop <- fread(input_file, header = F)
colnames(pop) <- c("#IID", "FID", "POP", "LAT", "LONG" )

# Sample ss_test individuals per deme
df_test <- pop %>% group_by(POP) %>% sample_n(ss_test)
test <- df_test[,1:2]

# Use the rest of the individuals for the GWAS panel
df_gwas <- anti_join(pop, df_test)
df_joint <- df_gwas[,1:2]


# Spilt GWAS panel based on Latitude
df_north <- df_gwas %>% filter(LAT >= 3)
gwas_north <- df_north[,1:2]
df_south <- df_gwas %>% filter(LAT < 3)
gwas_south <- df_south[,1:2]

# Write panel IDs to file
write.table(gwas_north, gwas_north_file, quote = F, col.names = F, row.names = F)
write.table(gwas_south, gwas_south_file, quote = F, col.names = F, row.names = F)
write.table(df_joint, gwas_joint_file, quote = F, col.names = F, row.names = F)
write.table(test, test_file, quote = F, col.names = F, row.names = F)

