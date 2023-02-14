# Format Tm into covariate file that will work with Plink 2 - this script includes 2 covariates (Tm plus popID/Latitude)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 7){stop("Usage: <pop file> <Tm> <fam file of GWAS> <test type> <output file> <PC file> <pc number>") }

# Parse args
pop_file = args[1]
Tm_file = args[2]
fam_file = args[3]
test_type = args[4]
output_file = args[5]
pc_file = args[6]

# Read in files
pops <- fread(pop_file, header = F)
fam <- fread(fam_file)
colnames(pops) <- c("#FID", "IID", "POP", "LAT", "LONG")

# Read in Tm
Tm <- fread(Tm_file)
Tm <- Tm$Tm

## Read in PCs
vecs <- fread(pc_file)
colnames(vecs)[1] <- "FID"
pcs <-  vecs[,3:ncol(vecs)]


if (test_type == "LAT") {

  # Format datafile
  tmp <- dplyr::inner_join(pops, fam) %>% select(c("#FID", "IID", "LAT"))
  df <- as.data.frame(cbind(tmp, Tm))
  df <- df[,c(1,2,4,3)]
  colnames(df) <- c("FID","IID", "Tm", "PopID")

} else if (test_type == "PS") {

  # Format datafile
  tmp <- dplyr::inner_join(pops, fam) %>% select(c("#FID", "IID", "POP"))
  df <- as.data.frame(cbind(tmp, Tm))
  df <- df[,c(1,2,4,3)]
  df <- df %>% mutate(PopID = case_when((POP == 25) ~ 1, (POP != 25) ~ 0)) %>% mutate(PopID = PopID - mean(PopID)) %>% select(c("#FID",	"IID",	"Tm", "PopID"))
  colnames(df) <- c("FID","IID", "Tm", "PopID")

} else {

  stop("Please enter acceptable test type: LAT, PS")

}

df <- cbind(df,pcs)

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
