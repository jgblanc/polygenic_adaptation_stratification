args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript simphenotype_ge_3.R <frequency file> <output_file> <seed>")}

suppressWarnings(suppressMessages({
library(data.table)
library(dplyr)
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
print(args[5])
set.seed(as.numeric(args[5]))
print(paste("The seed is", as.numeric(args[5])))

# load variant frequency file
p = fread(freq_file)
colnames(p)=c("chr","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("chr","ID","ALT_FREQS")]
p[, c("CHROM", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","ALT_FREQS")]
p$position = as.numeric(p$position)

# Function to randomly sample one variant per chromosome
sample.variant <- function(df) {
  return(sample_n(df,1))
}

# Sample one "causal" variant per chromosome
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

#let's calculate sigma2_g to confirm that the total genetic variance is indeed h2
sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, causal.variants$beta, causal.variants$ALT_FREQS))
print(paste0("The heritability is  ", sigma2_g))

#save the effect sizes to file and use plink2 to generate PRS
fwrite(causal.variants%>%
         mutate(ALT = "T")%>%
         select(ID,ALT,beta),
           effects_file,
       row.names=F,col.names=F,quote=F,sep="\t")


