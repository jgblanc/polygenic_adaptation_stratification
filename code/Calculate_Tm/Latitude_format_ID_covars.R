# Format Tm into covariate file that will work with Plink 2 - this script includes 2 covariates (Tm plus popID)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <pop file> <Tm> <fam file of GWAS> <output file> ") }

# Parse args
pop_file = args[1]
Tm_file = args[2]
fam_file = args[3]
output_file = args[4]

# Read in files
pops <- fread(pop_file, header = F)
fam <- fread(fam_file)
colnames(pops) <- c("#FID", "IID", "POP", "LAT", "LONG")

# Read in Tm
Tm <- fread(Tm_file)
Tm <- Tm$V1

# Format datafile
tmp <- dplyr::inner_join(pops, fam) %>% select(c("#FID", "IID", "LAT"))
df <- as.data.frame(cbind(tmp, Tm))
df <- df[,c(1,2,4,3)]
colnames(df) <- c("FID","IID", "Tm", "PopID")

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
