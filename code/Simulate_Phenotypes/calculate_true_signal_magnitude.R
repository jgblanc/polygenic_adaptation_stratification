args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simphenotype_ge_3.R <frequency file> <output_file> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(pgenlibr)
}))

# Popfile
pop_file=args[1]
print(paste("The pop file is", pop_file))

# Genotype prefix
geno_prefix = args[3]
print(paste("The test panel plink file prefix is", geno_prefix))

# Effect size file
es_file = args[2]
print(paste("The true effect size file is", es_file))

# Output file
out_file = args[4]
print(paste("The output file is", out_file))


# Function to read in genotype matrix for a set of variants
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat <- betas$BETA
  Z <- X %*% Bhat

  out <- cbind(Z)
  colnames(out) <- c("GV")

  return(out)
}

# Load effect sizes
betas <- fread(es_file)
colnames(betas) <- c("ID", "A1", "BETA")

# Load Genotypes
X <- read_genos(geno_prefix, betas)

# Calc GV
sscore <- pgs(X, betas)

# Read in fam file and add genetic values
fam <- fread(paste0(geno_prefix, ".psam"))
fam$GV <- sscore[,1]

# Read in popfile
pops <- fread(pop_file, header = F)
colnames(pops) <- c("#FID", "IID", "Pop")

# Compute average PGS in two test pops
df <- inner_join(pops, fam)
avgGV <- df %>% group_by(Pop) %>% summarise(avg_GV = mean(GV))
magnitude <- avgGV[1,2] - avgGV[2,2]
colnames(magnitude) <- "ts_magnitude"

# Write to output file
fwrite(magnitude, out_file,row.names=F,quote=F,sep="\t", col.names = T)

