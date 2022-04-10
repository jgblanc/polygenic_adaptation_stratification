args=commandArgs(TRUE)

if(length(args)<2){stop("<ids> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

id_file = args[1]
outfile = args[2]

ids <- fread(id_file)
test_pops <- unique(ids$Pop)
ids <- ids %>% mutate(Tvec = case_when(Pop == test_pops[1] ~ 1, Pop == test_pops[2] ~ 0))
ids <- ids %>% mutate(Tvec = Tvec - mean(Tvec))

fwrite(ids, outfile, row.names=F,quote=F,sep="\t", col.names = T)
