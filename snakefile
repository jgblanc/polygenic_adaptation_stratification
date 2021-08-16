CHR =[]
for i in range(0, 200):
  CHR.append(str(i))
CONFIG=["C1","C2"]
MODEL=["4PopSplit"]
REP = []
for i in range(1, 101):
  REP.append("F"+str(i))
HERITABILITY = ["seed-0"]
ENV = ["env-0.0","env-0.005", "env-0.01", "env-0.02", "env-0.03", "env-0.04", "env-0.06"]
#ENV = ["env-0.0"]
SIZE=2000
NUM_RESAMPLE=1000
PVALUE_THRESHOLD=1

def get_params(x):
  out = x.split("-")[1]
  return out

def get_seed(rep, h2, env):
  out1 = list(rep)[1]
  out2 = h2.split("-")[1]
  tmp = env.split("-")[1].split(".")[1]
  tmp_list = [int(i) for i in tmp]
  out3 = sum(tmp_list)
  return out1 + out2 + str(out3)

def get_seed1(rep, h2):
  out1 = list(rep)[1]
  out2 = h2.split("-")[1]
  return out1 + out2

rule all:
    input:
        expand("output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/Qx.txt",rep=REP, model = MODEL, h2 = HERITABILITY, env=ENV, config=CONFIG)

# Simluate Genotypes

rule simulate_genotypes_4popsplit:
    output:
        expand("output/Simulate_Genotypes/4PopSplit/{{rep}}/genos_{chr}.vcf", chr=CHR),
	      "output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    params:
      chr_num = len(CHR)
    shell:
        "python code/Simulate_Genotypes/generate_genotypes_4PopSplit.py \
	       --outpre output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
	       --chr {params.chr_num} \
       	       --Nanc 10000 \
	       --NA 10000 \
	       --NB 10000 \
	       --NC 10000 \
	       --ND 10000 \
  	     -a 2000 \
	       -b 2000 \
	       -c 2000 \
	       -d 2000 \
         -s1 22000 \
          -s2 11000 \
          -L 100000"

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
        """
        bcftools concat {input} -o output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/temp.vcf.gz -O z
        bcftools annotate --rename-chrs code/Simulate_Genotypes/convert_chr.txt output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/temp.vcf.gz -o {output} -O z
        rm output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/temp.vcf.gz
        """


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
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.psam",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pgen",
	"output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pvar",
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
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test-big \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/genos \
	      --rm-dup exclude-all
	            """

rule downsample_test:
    input:
      psam="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.psam",
          pgen="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pgen",
	      pvar="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pvar"
    output:
      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.psam",
      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pvar",
      "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test.pgen",
      "output/Simulate_Genotypes/{model}/{rep}/{config}/downsample.id"
    params:
      size = SIZE
    shell:
        """
	set +o pipefail;
        awk 'NR > 1' {input.psam} | cut -f 1,2 | sort -R | head -n {params.size} > output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/downsample.id

        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test-big \
        --keep output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/downsample.id \
        --make-pgen \
        --out output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test
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
	      psam=expand("output/Simulate_Genotypes/{model}/{rep}/genos.psam", model=MODEL, rep=REP, config=CONFIG),
	      big_psam=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.psam", model=MODEL, rep=REP, config=CONFIG),
	      big_pvar=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pvar", model=MODEL, rep=REP, config=CONFIG),
	      big_pgen=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test-big.pgen", model=MODEL, rep=REP, config=CONFIG),
        id=expand("output/Simulate_Genotypes/{model}/{rep}/{config}/downsample.id", model=MODEL, rep=REP, config=CONFIG)
    output:
        expand("output/Simulate_Genotypes/{model}/{rep}/ff.txt", model=MODEL, rep=REP)
    shell:
        """
	 touch {output}
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
	 rm {input.big_psam}
	 rm {input.big_pvar}
	 rm {input.big_pgen}
	 rm {input.id}
	 """

# Simluate Phenotypes

rule draw_effect_sizes:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.effects.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        seed = lambda wildcards: get_seed1(wildcards.rep, wildcards.h2)
    shell:
        "Rscript code/Simulate_Phenotypes/simgeffects.R {input} {output} {params.her} 0.4 {params.seed}"

rule generate_genetic_values:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pgen",
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.effects.txt"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.gvalue.sscore"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	            --out output/Simulate_Phenotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/genos-gwas_common.gvalue \
		          --score output/Simulate_Phenotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/genos-gwas_common.effects.txt cols=dosagesum,scoresums"

rule simulate_phenotype_4PopSplit:
    input:
        gvalues="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/genos-gwas_common.gvalue.sscore",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        "output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.phenos.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        en = lambda wildcards: get_params(wildcards.env),
        seed = lambda wildcards: get_seed(wildcards.rep,wildcards.h2,wildcards.env)
    shell:
        "Rscript code/Simulate_Phenotypes/simulate_phenotypes_4PopSplit_meanshift.R {input.gvalues} {input.pops} {output} {params.her} {params.en} {params.seed}"


# Run GWAS

rule gwas_no_correction:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.phenos.txt"
    output:
        "output/Run_GWAS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm allow-no-covars\
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common"



# Ascertain SNPs to include in PRS

rule pick_SNPS:
    input:
        causal_effect="output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    output:
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.betas",
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.p.betas",
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump_strat_only.R {input.causal_effect} output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common {params.pt} output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common"

# Get allele freq in test panel

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

# Calculate true breeding value

rule calc_true_gv:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        causal_effect="output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.effects.txt",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq"
    output:
        "output/PRS/{model}/{rep}/{config}/{h2}/genos-test_common.true.sscore",
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.causal_effect} cols=dosagesum,scoresums \
        --out output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/genos-test_common.true \
        """

## Include Tm as a covariate

# Generate Test Vector

rule make_test_vector:
    input:
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        fam="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam"
    output:
        "output/Calculate_Tm/{model}/{rep}/{config}/Tvec.txt"
    shell:
        "Rscript code/Calculate_Tm/make_tvec.R {input.pops} {input.fam} {output}"

# Project T using Plink2

rule proj_T:
    input:
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pgen",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.psam",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pgen",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        "output/Calculate_Tm/{model}/{rep}/{config}/Tvec.txt"
    params:
        n_minus_1 = int(SIZE)-1,
        col_start = 6,
        col_end = int(SIZE) + 4
    output:
        "output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenvec",
        "output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenval",
        "output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenvec.allele",
        "output/Calculate_Tm/{model}/{rep}/{config}/projection.sscore"
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-test_common \
       --pca allele-wts {params.n_minus_1} \
       --out output/Calculate_Tm/{wildcards.model}/{wildcards.rep}/{wildcards.config}/pca

        plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
       --score output/Calculate_Tm/{wildcards.model}/{wildcards.rep}/{wildcards.config}/pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
       --score-col-nums {params.col_start}-{params.col_end} \
       --out output/Calculate_Tm/{wildcards.model}/{wildcards.rep}/{wildcards.config}/projection
        """

# Calculate Tm using plink output

rule calc_Tm:
    input:
        vecs="output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenvec",
        vals="output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenval",
        proj="output/Calculate_Tm/{model}/{rep}/{config}/projection.sscore",
        tvec="output/Calculate_Tm/{model}/{rep}/{config}/Tvec.txt",
	allele="output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenvec.allele"
    output:
        "output/Calculate_Tm/{model}/{rep}/{config}/Tm.txt"
    shell:
        """
	Rscript code/Calculate_Tm/calc_Tm.R {input.vecs} {input.vals} {input.proj} {input.tvec} {output}

	rm {input.allele}
	"""

# Format Covariate file

rule format_covars:
    input:
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        fam="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        Tm="output/Calculate_Tm/{model}/{rep}/{config}/Tm.txt"
    output:
        "output/Calculate_Tm/{model}/{rep}/{config}/Tm_covars.txt"
    shell:
        "Rscript code/Calculate_Tm/format_covar.R {input.pops} {input.Tm} {input.fam} {output}"


# Re-run GWAS

rule gwas_Tm:
    input:
        genos="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.phenos.txt",
        Tm="output/Calculate_Tm/{model}/{rep}/{config}/Tm_covars.txt"
    output:
        "output/Run_GWAS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/{wildcards.model}/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm hide-covar \
        --covar {input.Tm} \
        --covar-col-nums 3 \
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common-Tm"

# Remake PRS

rule pick_SNPS_Tm:
    input:
        causal_effect="output/Simulate_Phenotypes/{model}/{rep}/{config}/{h2}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.pheno_strat.glm.linear"
    output:
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.betas",
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.p.betas",
        "output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump_strat_only.R {input.causal_effect} output/Run_GWAS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common-Tm {params.pt} output/PRS/{wildcards.model}/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.env}/genos-gwas_common-Tm"


# Do polygenic adaptation test

rule calc_lambdaT:
    input:
        vecs="output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenvec",
        vals="output/Calculate_Tm/{model}/{rep}/{config}/pca.eigenval",
        tvec="output/Calculate_Tm/{model}/{rep}/{config}/Tvec.txt",
    output:
        "output/Calculate_Tm/{model}/{rep}/{config}/Lambda_T.txt"
    shell:
        """
	Rscript code/Calculate_Tm/calc_lambdaT.R {input.vecs} {input.vals} {input.tvec} {output}
	"""

rule calc_Va:
    input:
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq",
        c="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.betas",
        cp="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.p.betas",
        nc="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common.nc.betas"
    output:
        "output/PGA_test/{model}/{rep}/{config}/{h2}/{env}/Va.txt"
    shell:
      """
          Rscript code/PGA_test/calc_Va_strat.R {input.freq} {input.c} {input.cp} {input.nc} {output}
	      """

rule calc_Va_Tm:
    input:
        freq="output/Simulate_Genotypes/{model}/{rep}/{config}/genos-test_common.afreq",
        c="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.betas",
        cp="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.p.betas",
        nc="output/PRS/{model}/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.nc.betas"
    output:
        "output/PGA_test/{model}/{rep}/{config}/{h2}/{env}/Va-Tm.txt"
    shell:
      """
          Rscript code/PGA_test/calc_Va_strat.R {input.freq} {input.c} {input.cp} {input.nc} {output}
	      """

rule Calc_Qx:
    input:
        c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.betas",
        cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common.c.p.betas",
        nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common.nc.betas",
        c_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.betas",
        cp_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.c.p.betas",
        nc_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{env}/genos-gwas_common-Tm.nc.betas",
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        lambda_T="output/Calculate_Tm/4PopSplit/{rep}/{config}/Lambda_T.txt",
        Va="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/Va.txt",
        Va_Tm="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/Va-Tm.txt",
        true="output/PRS/4PopSplit/{rep}/{config}/{h2}/genos-test_common.true.sscore",
        Tvec="output/Calculate_Tm/4PopSplit/{rep}/{config}/Tvec.txt",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        qx="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/Qx.txt",
        pgs="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/PGS.txt"
    params:
        num=NUM_RESAMPLE
    shell:
      """
          Rscript code/PGA_test/calc_Qx_4PopSplit.R {input.c} {input.cp} {input.nc} {input.c_Tm} {input.cp_Tm} {input.nc_Tm} output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common {input.lambda_T} {input.Va} {input.Va_Tm} {input.true} {input.Tvec} {input.pops} {params.num} {output.qx} {output.pgs}
	      """

# Remove non-end files

#rule delete_files:
#    input:
#        qx=expand("output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/Qx.txt", model=MODEL, rep=REP, config=CONFIG, #h2=HERITABILITY, env=ENV),
#        pgs=expand("output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{env}/PGS.txt", model=MODEL, rep=REP, config=CONFIG, h2=HERITABILITY, env=ENV)
#    output:
#        "output/PGA_test/4PopSplit/{rep}/finish.txt"
#    shell:
#        """
#	 touch {output}
#	 rm -r output/Simulate_Genotypes/4PopSplit/{wildcards.rep}
#	 rm -r output/Simulate_Phenotypes/4PopSplit/{wildcards.rep}
#	 rm -r output/Calculate_Tm/4PopSplit/{wildcards.rep}
#	 rm -r output/Run_GWAS/4PopSplit/{wildcards.rep}
#	 rm -r output/PRS/4PopSplit/{wildcards.rep}
#	 """

