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
colnames(fam) <- c("FID", "IID", "SEX")

# Make mean centered Test vec (PopID)
pops <- fread(pop_file, header = T)
pop <- dplyr::inner_join(pops, fam) %>% dplyr::select(c("FID", "IID", "POP"))
test_pops <- unique(pop$POP)
pop <- pop %>% mutate(Tvec = case_when(POP == test_pops[1] ~ 1, POP == test_pops[2] ~ 0)) %>% mutate(Tvec = Tvec - mean(Tvec))

# Write to file
write.table(pop, output_file,row.names=F,quote=F,sep="\t", col.names = T)
