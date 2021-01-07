
rule simulate_genotypes_4popsplit:
    output:
        "output/Simulate_Genotypes/4PopSplit/test.vcf",
	"output/Simulate_Genotypes/4PopSplit/test.pop"
    shell:
        "python code/Simulate_Genotypes/generate_genotypes_4PopSplit.py \
	--outpre output/Simulate_Genotypes/4PopSplit/test"

rule format_VCF:
    input:
        "output/Simulate_Genotypes/4PopSplit/test.vcf"
    output:
        "output/Simulate_Genotypes/4PopSplit/test.ids.vcf.gz"
    shell:
        """
	head -n6 {input} > output/Simulate_Genotypes/4PopSplit/header.txt,
	cat output/Simulate_Genotypes/4PopSplit/header.txt <(cat output/Simulate_Genotypes/4PopSplit/test.vcf | awk -v OFS="\t" 'NR>6 {{$3=$1"_"$2"_A_T";$4="A"; $5="T"; print ;}}') | bgzip > {output}
	bcftools index {output}
	"""

rule LD_clumping_GWAS_atlas:
    input:
        "output/GWAS_ATLAS/parsed_gwas/{trait}_parsed.txt"
    output:
        "output/GWAS_ATLAS/clumped/{trait}.clumped"
    shell:
        "code/plink \
        --noweb \
        --bfile data/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01_VARID \
        --clump {input} \
        --clump-field P \
        --clump-p1 1 \
        --clump-p2 1 \
        --clump-r2 0.5 \
        --clump-kb 250 \
        --out output/GWAS_ATLAS/clumped/{wildcards.trait}"

rule Extract_clumped_SNPs_GWAS_ATLAS:
    input:
        "output/GWAS_ATLAS/clumped/{trait}.clumped"
    output:
        "output/GWAS_ATLAS/clumped/{trait}_SNPs.txt"
    shell:
        "awk '{{ print $3}}' {input} > {output}"

rule LD_pruning_GWAS_atlas:
    input:
        "output/GWAS_ATLAS/parsed_gwas/{trait}_{threshold}_parsed.txt"
    output:
        "output/GWAS_ATLAS/pruned/{trait}_{threshold}.prune.in"
    shell:
        """
	cut -f 1 -d' ' {input} > all_ss.snps
	code/plink \
    	--noweb \
    	--bfile data/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01/1000G_20101123_v3_GIANT_chr1_23_minimacnamesifnotRS_CEU_MAF0.01_VARID \
   	--extract all_ss.snps \
    	--make-bed \
    	--out all_ss_plink
	code/plink \
    	--bfile all_ss_plink \
    	--indep-pairwise 50 5 0.05 \
    	--noweb \
    	--out output/GWAS_ATLAS/pruned/{wildcards.trait}_{wildcards.threshold}
	rm all_ss*
	"""

rule LD_pruning_STRAT:
    input:
        "data/STRAT/chr1_EUR_{MAF}.eigenvec.var.DA.txt"
    output:
        "output/STRAT/pruned/chr1_EUR_{MAF}.eigenvec.var.DA.prune.in"
    shell:
        """
        cut -f2,25 -d',' {input} > all_ss.temp
	awk '$2 != "NA"' FS=',' all_ss.temp | cut -f1 -d',' > all_ss.snps #Pick only D/A SNPs and get rsID 
        code/plink \
        --noweb \
        --bfile data/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR \
        --extract all_ss.snps \
        --make-bed \
        --out all_ss_plink
        code/plink \
        --bfile all_ss_plink \
        --indep-pairwise 50 5 0.05 \
        --noweb \
        --out output/STRAT/pruned/chr1_EUR_{wildcards.MAF}.eigenvec.var.DA
	rm all_ss*
        """