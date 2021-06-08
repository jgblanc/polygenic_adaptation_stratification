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

# Make standardized test vec
pops <- fread(pop_file, header = F)
colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))
Tvec <- pop$Lat
Tvec <- (Tvec-mean(Tvec))
#Tvec <- ctvec/sqrt(sum(ctvec^2)/length(ctvec))

# Write to file
write.table(Tvec, output_file,row.names=F,quote=F,sep="\t", col.names = F)
