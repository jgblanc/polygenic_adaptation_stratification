# Get list of SNPs that are at 5% MAF or greater in both the Test and GWAS panel
library(data.table)

# GWAS .freq files
north <- fread(snakemake@input[[1]], header = T)
south <- fread(snakemake@input[[2]], header = T)

# Test panel .freq
test <- fread(snakemake@input[[3]], header = T)

# Subset variants
n <- subset(north, north$ALT_FREQS > 0.05 & north$ALT_FREQS < 0.95)
s <- subset(south, south$ALT_FREQS > 0.05 & south$ALT_FREQS < 0.95)
t <- subset(test, test$ALT_FREQS > 0.05 & test$ALT_FREQS < 0.95)

# Merge to include only variants passing frequency filter in both panels
dat <- merge(n,s, by=c("#CHROM","ID"))
dat2 <- merge(dat,t, by=c("#CHROM","ID"))

write.table(dat$ID, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
