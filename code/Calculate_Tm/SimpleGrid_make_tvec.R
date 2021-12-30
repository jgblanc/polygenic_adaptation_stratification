# Make Test Vector

library(dplyr)
library(data.table)

args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <pop file> <fam file> <type of test> <output file> ") }

# Parse args
pop_file = args[1]
fam_file = args[2]
test_type = args[3]
output_file = args[4]

# Read in Fam file
fam <- fread(fam_file)
colnames(fam) <- c("IID", "FID", "SEX")

if (test_type == "LAT") {

  print(test_type)

  pops <- fread(pop_file, header = F)
  colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
  pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))

  # Test vector is latitude
  Tvec <- pop$Lat

  # Mean center
  Tvec <- (Tvec-mean(Tvec))

} else if (test_type == "PS") {

  print(test_type)

  pops <- fread(pop_file, header = F)
  colnames(pops) <- c("IID", "FID", "Pop", "Lat", "Long")
  pop <- dplyr::inner_join(pops, fam, by = c("IID"= "IID"))

  # Test vector is deme 25 vs everyone else (mean center)
  pop <- pop %>% mutate(T1 = case_when(Pop == 25 ~ 1, Pop != 25 ~ 0)) %>% mutate(T1 = T1 - mean(T1))
  Tvec <- pop$T1

} else {
  stop("Please enter acceptable test type: LAT, PS")
}


# Write to file
write.table(Tvec, output_file,row.names=F,quote=F,sep="\t", col.names = F)
