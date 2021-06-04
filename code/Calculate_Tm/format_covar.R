# Format Tm into covariate file that will work with Plink 2

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

# Read in Tm
Tm <- fread(Tm_file)
Tm <- Tm$V1

# Format datafile
tmp <- dplyr::inner_join(pops, fam, by = c("V1"= "IID")) %>% select("V2", "V3")
df <- as.data.frame(cbind(tmp$V2, tmp$V2, Tm))
colnames(df) <- c("FID","IID", "Tm")

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
