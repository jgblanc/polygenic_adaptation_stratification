# Format Tm into covariate file that will work with Plink 2 - this script includes 2 covariates (Tm plus popID)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <pop file> <Tm> <fam file of GWAS> <output file> <eigenvec>") }

# Parse args
pop_file = args[1]
Tm_file = args[2]
fam_file = args[3]
output_file = args[4]
#evec_file = args[5]

# Read in files
pops <- fread(pop_file, header = T)
fam <- fread(fam_file)

# Read in Tm
Tm <- fread(Tm_file)

# Format datafile
tmp <- dplyr::inner_join(pops, fam, by = c("IID"= "IID")) %>% select("IID", "POP")
df <- as.data.frame(cbind(tmp$IID, tmp$IID, Tm, tmp$POP))
colnames(df) <- c("FID","IID", "Tm", "PopID")
pops <- unique(df$PopID)
df <- df %>% mutate(PopID = case_when((PopID == pops[1]) ~ 1, (PopID == pops[2]) ~ 0))

## Read in eigenvec files
#vecs <- fread(evec_file)
#colnames(vecs)[1] <- "FID"
#df <- inner_join(df, vecs)
df$PC1 <- rnorm(nrow(df))
print(head(df))

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
