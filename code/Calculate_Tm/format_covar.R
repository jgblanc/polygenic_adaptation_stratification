# Format Tm into covariate file that will work with Plink 2

library(dplyr)
library(data.table)
library(rhdf5)

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

# Make Tm
#file <- h5file(Tm_file)
#file$ls(recursive = T)
#n <- file[['output']]$maxdims
#Tm <- rep(0, n)
#for (i in 1:n) {
#  Tm[i] = file[['output']][i]
#}
Tm <- h5read(Tm_file, "output/")

# Format datafile
tmp <- dplyr::inner_join(pops, fam, by = c("V1"= "IID")) %>% select("V2", "V3")
df <- as.data.frame(cbind(tmp$V2, tmp$V2, Tm))
colnames(df) <- c("FID","IID", "Tm")

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
