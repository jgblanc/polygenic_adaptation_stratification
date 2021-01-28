CHR=["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
#CHR=["0", "1"]
CONFIG=["C1","C2"]
MODEL=["4PopSplit"]
REP=["V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"]


rule all:
    input:
        expand("output/PRS/{model}/{rep}/{config}/genos-test_common.c.p.sscore", model=MODEL, rep=REP, config=CONFIG),
        expand("output/PRS/{model}/{rep}/{config}/genos-test_common.true.sscore", model=MODEL, rep=REP, config=CONFIG)

# Simluate Genotypes

rule simulate_genotypes_4popsplit:
    output:
        expand("output/Simulate_Genotypes/4PopSplit/{{rep}}/genos_{chr}.vcf", chr=CHR),
	      "output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    shell:
        "python code/Simulate_Genotypes/generate_genotypes_4PopSplit.py \
	       --outpre output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
	       --chr 20 \
	       --Nanc 40000 \
	       -a 5000 \
	       -b 5000 \
	       -c 5000 \
	       -d 5000"

rule format_VCF:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/genos_{chr}.vcf"
    output:
        gz="output/Simulate_Genotypes/{model}/{rep}/genos_{chr}.ids.vcf.gz"
	      #csi="output/Simulate_Genotypes/{model}/{rep}/genos_{chr}.ids.vcf.gz.csi"
    shell:
        """
	      head -n6 {input} > output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/header_{wildcards.chr}.txt
	      cat output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/header_{wildcards.chr}.txt <(cat output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/genos_{wildcards.chr}.vcf | awk -v OFS="\t" 'NR>6 {{$3=$1"_"$2"_A_T";$4="A"; $5="T"; print ;}}') | bgzip > {output.gz}
	      #bcftools index {output.gz}
	      rm output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/header_{wildcards.chr}.txt
	      """

rule concat_vcfs:
    input:
        expand("output/Simulate_Genotypes/{{model}}/{{rep}}/genos_{chr}.ids.vcf.gz", chr=CHR)
    output:
        "output/Simulate_Genotypes/{model}/{rep}/genos.ids.vcf.gz"
    shell:
        "bcftools concat {input} -o {output} -O z"

rule convert_vcf_to_plink:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/genos.ids.vcf.gz"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/genos.psam",
	"output/Simulate_Genotypes/{model}/{rep}/genos.pgen",
      	"output/Simulate_Genotypes/{model}/{rep}/genos.pvar"
    shell:
        "plink2 \
        --double-id \
        --make-pgen \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/genos \
        --vcf {input}"

rule create_panels_4PopSplit:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/genos.pop"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/C1/ids.gwas",
	      "output/Simulate_Genotypes/{model}/{rep}/C1/ids.test",
	      "output/Simulate_Genotypes/{model}/{rep}/C2/ids.gwas",
	      "output/Simulate_Genotypes/{model}/{rep}/C2/ids.test"
    script:
        "code/Simulate_Genotypes/split_gwas-test_4PopSplit.R"

rule split_into_test_gwas:
    input:
        gwas="output/Simulate_Genotypes/{model}/{rep}/{config}/ids.gwas",
	      test="output/Simulate_Genotypes/{model}/{rep}/{config}/ids.test",
	      psam="output/Simulate_Genotypes/{model}/{rep}/genos.psam",
	      pvar="output/Simulate_Genotypes/{model}/{rep}/genos.pvar",
	      pgen="output/Simulate_Genotypes/{model}/{rep}/genos.pgen"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.psam",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pgen",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pvar",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.psam",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pgen",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pvar"
    shell:
        """
	      plink2 \
	      --keep {input.gwas} \
	      --make-pgen \
	      --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/genos \
	      --rm-dup exclude-all

        plink2 \
        --keep {input.test} \
        --make-pgen \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/genos \
	      --rm-dup exclude-all
	      """

rule get_variant_freq:
    input:
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.psam",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pvar",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pgen"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.afreq",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.afreq"
    shell:
        """
        plink2 \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test \
	      --freq \
	      --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test

	      plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas \
        --freq \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas
        """

rule get_common_snp_list:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.afreq",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.afreq"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/common_snp_ids.txt"
    script:
        "code/Simulate_Genotypes/get_common_snp_list.R"

rule remake_panels_with_common_snps:
    input:
        common_id="output/Simulate_Genotypes/{model}/{rep}/{config}/common_snp_ids.txt",
	      test_psam="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.psam",
	      test_pvar="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pvar",
	      test_pgen="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pgen",
	      gwas_psam="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.psam",
	      gwas_pvar="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pvar",
	      gwas_pgen="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pgen"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pvar",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pgen",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pvar",
	      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pgen"
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
	      --make-pgen

        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	      --make-pgen
        """

rule common_snp_freq:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pgen"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	      --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	      --freq"

rule aggregate_genotypes:
    input:
        frq=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq", model=MODEL, rep=REP, config=CONFIG),
	genos=expand("output/Simulate_Genotypes/{model}/{rep}/genos_{chr}.vcf", chr=CHR, rep=REP, config=CONFIG, model=MODEL),
	gz_chr=expand("output/Simulate_Genotypes/{model}/{rep}/genos_{chr}.ids.vcf.gz", chr=CHR, model=MODEL, rep=REP),
	frq_test=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.afreq", model=MODEL, rep=REP, config=CONFIG),
	frq_gwas=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.afreq", model=MODEL, rep=REP, config=CONFIG),
	gz=expand("output/Simulate_Genotypes/{model}/{rep}/genos.ids.vcf.gz", model=MODEL, rep=REP),
	gwas_pgen=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pgen", model=MODEL, rep=REP, config=CONFIG),
	gwas_pvar=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.pvar", model=MODEL, rep=REP, config=CONFIG),
	gwas_psam=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas.psam", model=MODEL, rep=REP, config=CONFIG),
	test_pgen=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pgen", model=MODEL, rep=REP, config=CONFIG),
	test_pvar=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pvar", model=MODEL, rep=REP, config=CONFIG),
	test_psam=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.psam", model=MODEL, rep=REP, config=CONFIG),
	pgen=expand("output/Simulate_Genotypes/{model}/{rep}/genos.pgen", model=MODEL, rep=REP, config=CONFIG),
	pvar=expand("output/Simulate_Genotypes/{model}/{rep}/genos.pvar", model=MODEL, rep=REP, config=CONFIG),
	psam=expand("output/Simulate_Genotypes/{model}/{rep}/genos.psam", model=MODEL, rep=REP, config=CONFIG) 
    shell:
        """
	echo {input.frq}
	rm {input.genos}
	rm {input.frq_test}
	rm {input.frq_gwas}
	rm {input.gz}
	rm {input.gz_chr}
	rm {input.gwas_pgen}
	rm {input.gwas_pvar}
	rm {input.gwas_psam}
	rm {input.test_pgen}
	rm {input.test_pvar}
	rm {input.test_psam}
	rm {input.pgen}
	rm {input.pvar}
	rm {input.psam}
	"""

# Simluate Phenotypes

rule draw_effect_sizes:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.effects.txt"
    shell:
        "Rscript code/Simulate_Phenotypes/simgeffects.R {input} {output} 0.8 0.4 12"

rule generate_genetic_values:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pgen",
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.effects.txt"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.gvalue.sscore"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	      --out output/Simulate_Phenotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common.gvalue \
	      --score output/Simulate_Phenotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common.effects.txt cols=dosagesum,scoresums"

rule simulatate_phenotype_4PopSplit:
    input:
        gvalues="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/genos-gwas_common.gvalue.sscore",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.phenos.txt"
    shell:
        "Rscript code/Simulate_Phenotypes/simulate_phenotype_4PopSplit.R {input.gvalues} {input.pops} {output} 0.8 12"

rule aggregate_phenotypes:
    input:
        expand("output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.phenos.txt", model=MODEL, rep=REP, config=CONFIG)
    shell:
        "echo {input}"

# Run GWAS

rule gwas_no_correction:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.phenos.txt"
    output:
        "output/Run_GWAS/{model}/{rep}/{config}/genos-gwas_common.pheno_random.glm.linear",
        "output/Run_GWAS/{model}/{rep}/{config}/genos-gwas_common.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm allow-no-covars \
        --pheno {input.pheno} \
        --pheno-name pheno_random,pheno_strat \
        --out output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common"

rule aggregate_gwas:
    input:
        expand("output/Run_GWAS/{model}/{rep}/{config}/genos-gwas_common.pheno_random.glm.linear", model=MODEL, rep=REP, config=CONFIG)
    shell:
        "echo {input}"

# PRS

rule pick_SNPS:
    input:
        causal_effect="output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.effects.txt",
        gwas_random="output/Run_GWAS/{model}/{rep}/{config}/genos-gwas_common.pheno_random.glm.linear",
        gwas_strat="output/Run_GWAS/{model}/{rep}/{config}/genos-gwas_common.pheno_strat.glm.linear"
    output:
        "output/PRS/{model}/{rep}/{config}/genos-gwas_common.c.betas",
        "output/PRS/{model}/{rep}/{config}/genos-gwas_common.c.p.betas",
        "output/PRS/{model}/{rep}/{config}/genos-gwas_common.nc.betas"
    shell:
        "Rscript code/PRS/clump.R {input.causal_effect} output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common 5e-4 output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common"

rule test_snp_freq:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pgen"
    output:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
	      --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
	      --freq"

rule calc_prs:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        c="output/PRS/{model}/{rep}/{config}/genos-gwas_common.c.betas",
        cp="output/PRS/{model}/{rep}/{config}/genos-gwas_common.c.p.betas",
        nc="output/PRS/{model}/{rep}/{config}/genos-gwas_common.nc.betas",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq"
    output:
        "output/PRS/{model}/{rep}/{config}/genos-test_common.c.sscore",
        "output/PRS/{model}/{rep}/{config}/genos-test_common.c.p.sscore",
        "output/PRS/{model}/{rep}/{config}/genos-test_common.nc.sscore"
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.c} cols=dosagesum,scoresums \
        --out output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common.c \
        --score-col-nums 3,4

        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.cp} cols=dosagesum,scoresums \
        --out output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common.c.p \
        --score-col-nums 3,4

        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.nc} cols=dosagesum,scoresums \
        --out output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common.nc \
        --score-col-nums 3,4
        """

rule calc_true_gv:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        causal_effect="output/Simulate_Phenotypes/{model}/{rep}/{config}/genos-gwas_common.effects.txt",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq"
    output:
        "output/PRS/{model}/{rep}/{config}/genos-test_common.true.sscore",
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.causal_effect} cols=dosagesum,scoresums \
        --out output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common.true \
        """

rule aggregate_prs:
    input:
        expand("output/PRS/{model}/{rep}/{config}/genos-test_common.c.p.sscore", model=MODEL, rep=REP, config=CONFIG),
        expand("output/PRS/{model}/{rep}/{config}/genos-test_common.true.sscore", model=MODEL, rep=REP, config=CONFIG)
    shell:
        "echo {input}"

