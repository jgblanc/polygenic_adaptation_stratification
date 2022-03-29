args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript draw_effects_sizes.R <frequency file> <output_file> <heritability> <alpha> <test panel genotype prefix> <popfile> <probability of flipping effect size>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
  library(tidyr)
  library(Hmisc)
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

#genotypes prefix
geno_prefix = args[5]

# pop file
popfile = args[6]

# probability effect size is positive given pC - pD is positive
prob = as.numeric(args[7])
print(paste("Prob is", prob))

# load variant frequency file
p = fread(freq_file)
freq <- p %>% dplyr::select(ID, ALT, ALT_FREQS) %>% separate(ID, c("CHR", "POS"), ":")
freq <- cbind(freq, p$ID)
colnames(freq) <- c("CHROM", "POS", "ALT", "ALT_FREQS", "ID")


# Choose one variant every 300 sites to be causal
sample.variant <- function(df) {
  nsnp <- length(df$POS)
  reps <- round(nsnp / 300)
  df$group <- as.numeric(cut2(as.numeric(df$POS), g = reps))
  out <- df %>% group_by(group) %>% sample_n(1)
  out <- out[,c(1,2,3,4,5)]
  return(out)
}

#carry this out grouped by chromosome
df = freq %>% filter(CHROM == 1)
nchrms <- length(unique(freq$CHROM))
causal.variants <- sample.variant(as.data.table(df))
if(nchrms > 1) {
  for (i in 2:nchrms){
    df = freq %>% filter(CHROM == i)
    out <- sample.variant(as.data.table(df))
    causal.variants <- rbind(causal.variants, out)
  }
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

#let's calculate sigma2_g to confirm that the total genetic variance is indeed 0.3
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
pop <- fread(popfile, header = T)
pop <- pop %>% select(Sample, Population)
colnames(pop) <- c("#IID", "POP")

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
G[is.na(G)] <- 0

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
         select(ID,ALT,beta),
           effects_file,
       row.names=F,col.names=F,quote=F,sep="\t")


