args=commandArgs(TRUE)

if(length(args)<8){stop("Rscript draw_effects_sizes_4PopSplit.R <frequency file> <output_file> <heritability> <alpha> <test panel genotype prefix> <popfile> <probability of flipping effect size> <seed>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
  library(tidyr)
}))

#frequency file
freq_file=args[1]
print(paste("The frequency file is",freq_file))

#output file - genetic effects
effects_file = args[2]
print(paste("The output file is",effects_file))

#heritability
h2 = as.numeric(args[3])
print(paste("The heritability is",h2))

#alpha
alpha = as.numeric(args[4])
print(paste("Alpha is", alpha))

#random seed
seed = as.numeric(args[5])
set.seed(seed)

#genotypes prefix
geno_prefix = args[6]

# pop file
popfile = args[7]

# probability effect size is positive given pC - pD is positive
prob = as.numeric(args[8])
print(paste("Prob is", prob))

# load variant frequency file
p = fread(freq_file)
colnames(p)=c("chr","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("chr","ID","ALT_FREQS")]
p[, c("CHROM", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","ALT_FREQS")]
p$position = as.numeric(p$position)


sample.variant <- function(df) {
  return(sample_n(df,1))
}

#carry this out grouped by chromosome
df = p %>% filter(CHROM == 1)
nchrms <- length(unique(p$CHROM))
causal.variants <- sample.variant(as.data.table(df))
for (i in 2:nchrms){
  df = p %>% filter(CHROM == i)
  out <- sample.variant(as.data.table(df))
  causal.variants <- rbind(causal.variants, out)
}


#Now generate the effect sizes from these variants
#calculate the independent component of variance required
sigma2_l = h2 / sum( sapply( causal.variants$ALT_FREQS,function(x){
  beta= ( 2*x*(1-x)) ^ (1-alpha)
  return(beta)
}))

#sample maf-dependent effects using the model above
causal.variants$beta = sapply( causal.variants$ALT_FREQS , function(x){
  beta = rnorm( 1 , mean = 0, sd = sqrt(sigma2_l * (2*x*(1-x))^-alpha ))
})
causal.variants$beta <- abs(causal.variants$beta)

#let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.8
sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, causal.variants$beta, causal.variants$ALT_FREQS))
print(sigma2_g)

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

# Read in population ID info
pop <- fread(popfile, header = F)
colnames(pop) <- c("FID", "IID", "POP")

# Read in fam file
fam <- fread(paste0(geno_prefix, ".psam"))

# Get only test individuals
test_inds <- inner_join(fam, pop)

# Get number of individuals in each population
n1 <- as.numeric(count(test_inds, POP)[1,2])
n2 <- as.numeric(count(test_inds, POP)[2,2])
print(c(n1, n2))

# Read in genotype matrix for causal variants
G <- read_genos(geno_prefix, causal.variants[,"ID"])

# Calculate population specific allele frequeny
p1 <- colMeans(G[1:n1,])/2
p2 <- colMeans(G[(n1+1):(n1+n2),])/2

# Get allele frequency difference
diff <- p1 - p2
print(max(diff))

# Create correlation between effect size and pop ID
for (i in 1:nrow(causal.variants)){
  b <- causal.variants[i,"beta"]
  if (diff[i] >= 0) {
    causal.variants[i,"beta"] <- sample(c(-1, 1),1, prob = c((1-prob), prob)) * b
  } 
}

# Print probability of positive beta
print(sum(causal.variants$beta > 0)/length(causal.variants$beta))

#save the effect sizes to file and use plink2 to generate PRS
fwrite(causal.variants%>%
         mutate(ALT = "T")%>%
         select(ID,ALT,beta),
           effects_file,
       row.names=F,col.names=F,quote=F,sep="\t")


