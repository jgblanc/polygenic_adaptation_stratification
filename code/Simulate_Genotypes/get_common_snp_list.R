# Get list of SNPs that are preset in both the Test and GWAS panel
library(data.table)

tp <- fread(snakemake@input[[1]], header = T)
gp <- fread(snakemake@input[[2]], header = T)

t <- subset(tp, tp$ALT_FREQS > 0.05 & tp$ALT_FREQS < 0.95)
g <- subset(gp, gp$ALT_FREQS > 0.05 & gp$ALT_FREQS < 0.95)

dat <- merge(t,g, by=c("#CHROM","ID"))

write.table(dat$ID, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
