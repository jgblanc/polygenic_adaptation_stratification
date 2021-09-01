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
pop <- pop %>% mutate(T1 = case_when(Pop == 7 ~ 1, Pop != 7 ~ 0)) %>% mutate(T1 = T1 - mean(T1))
pop <- pop %>% mutate(T2 = case_when(Pop == 22 ~ 1, Pop != 22 ~ 0)) %>% mutate(T2 = T2 - mean(T2))
pop <- pop %>% mutate(T3 = case_when(Pop == 32 ~ 1, Pop != 32 ~ 0)) %>% mutate(T3 = T3 - mean(T3))


# Write to file
df <- pop %>% select(T1, T2, T3)
fwrite(df,output_file, col.names=T,row.names=F,quote=F,sep="\t")
