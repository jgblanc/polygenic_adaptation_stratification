# Split 1KG individuals into gwas and test panels

args=commandArgs(TRUE)

if(length(args)!=4){stop("Rscript split_gwas-test.R <sample info> <fam file> <outfile gwas> <outifle test>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

popfile = args[1] # Sample info
famfile = args[2]
gwas_outfile = args[3]
test_outfile = args[4]

# Read files
pop <- fread(popfile)
fam <- fread(famfile)

# Join based on sample ID
df <- inner_join(pop, fam, by = c("Sample"= "#IID")) %>% dplyr::select(Sample, Population)

# Pick 50 FIN and 50 IBS to be test panel
test <- df %>% filter(Population == "FIN" | Population == "IBS") %>% group_by(Population) %>% sample_n(50)

# Use the rest of individuals for GWAS panel
gwas <- anti_join(df, test)

# Save
colnames(test) <- c("#IID", "Pop")
fwrite(test, test_outfile, quote = F, col.names = T, row.names = F, sep = "\t")
colnames(gwas) <- c("#IID", "Pop")
fwrite(gwas, gwas_outfile, quote = F, col.names = T, row.names = F, sep = "\t")

