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

# Make mean centered Test vec (PopID)
pops <- fread(pop_file, header = F)
pop <- dplyr::inner_join(pops, fam, by = c("V1"= "IID")) %>% select("V2", "V3")
test_pops <- unique(pop$V3)
pop <- pop %>% mutate(tvec = case_when(V3 == test_pops[1] ~ 1, V3 == test_pops[2] ~ 0))
tvec <- pop$tvec
Tvec <- (tvec-mean(tvec))

# Add test vec to new column
colnames(fam) <- c("FID", "IID", "Tvec")
fam$Tvec <- Tvec

# Write to file
write.table(fam, output_file,row.names=F,quote=F,sep="\t", col.names = T)
