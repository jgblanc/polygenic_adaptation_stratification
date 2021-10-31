# Make Test Vector

library(dplyr)
library(data.table)

args=commandArgs(TRUE)

if( length(args) != 3){stop("Usage: <pop file> <fam file> <output file> ") }

# Parse args
pop_file = args[1]
fam_file = args[2]
output_file = args[3]

# Read in Fam file
fam <- fread(fam_file)
colnames(fam) <- c("IID", "FID", "SEX")

# Join files
pops <- fread(pop_file, header = F)
colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))

# Make Tvecs
pop <- pop %>% mutate(T1 = case_when(Pop == 25 ~ 1, Pop != 25 ~ 0)) %>% mutate(T1 = T1 - mean(T1))

# Write to file
df <- pop %>% select(T1)
fwrite(df,output_file, col.names=T,row.names=F,quote=F,sep="\t")
