# Format Tm into covariate file that will work with Plink 2
suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 6){stop("Usage: <pop file> <Tm> <fam file of GWAS> <test type> <output file> <PC file>") }

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

} else if (test_type == "DIAG") {

  id_diag <- c(0,7,14, 21, 28, 35)
  # Format datafile
  tmp <- dplyr::inner_join(pops, fam) %>% select(c("#FID", "IID", "POP"))
  df <- as.data.frame(cbind(tmp, Tm))
  df <- df[,c(1,2,4,3)]
  df <- df %>% mutate(PopID = case_when((POP %in% id_diag) ~ 1, !(POP %in% id_diag) ~ 0)) %>% mutate(PopID = PopID - mean(PopID)) %>% select(c("#FID",	"IID",	"Tm", "PopID"))
  colnames(df) <- c("FID","IID", "Tm", "PopID")


}else {

  stop("Please enter acceptable pheno type: LAT, PS, DIAG")

}

df <- cbind(df,pcs)

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
