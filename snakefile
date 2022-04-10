CHR =[]
for i in range(0, 200):
  CHR.append(str(i))
REP = []
for i in range(1, 101):
  REP.append("B"+str(i))
HERITABILITY = ["h2-0"]
PHENO = ["LAT"]
ENV = ["env-0.5"]
SS_TEST = 20 # Number indiviuals per deme in test panel
SIZE = SS_TEST * 36
GWAS_SIZE = 60 * 36
PVALUE_THRESHOLD = 1
NUM_RESAMPLE = 1000


wildcard_constraints:
    rep="[A-Z]\d+",
    config="C1",
    h2="h2-[0-1]",
    pheno="[A-Z]*",
    env="env-[0-9].[0-9]*",
    test="[A-Z]*"

def get_params(x):
  out = x.split("-")[1]
  return out

def get_seed(rep, h2, pheno, env):
  rep = str(sum(list(map(int,list(rep)[1::]))))
  h2 = h2.split("-")[1]
  pheno_list = list(pheno)
  pheno_list = [ord(i) for i in pheno_list]
  pheno = sum(pheno_list)
  env_list = env.split("-")[1].split(".")[1]
  env_list = [int(i) for i in env_list]
  env = sum(env_list)
  out = rep + h2 + str(pheno) + str(env)
  return out

def get_seed_msprime(rep):
  out = int(''.join(list(rep)[1::])) * 1000
  return out


rule all:
    input:
        expand("output/Run_GWAS/SimpleGrid/{rep}/South/{h2}/{pheno}/{env}/genos-gwas_common.pheno_strat.glm.linear", rep=REP, h2=HERITABILITY, env=ENV, pheno=PHENO)

# Simluate Genotypes

rule simulate_genotypes_SimpleGrid:
    output:
        expand("output/Simulate_Genotypes/SimpleGrid/{{rep}}/genos_{chr}.vcf", chr=CHR),
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pop"
    params:
      chr_num = len(CHR),
      seed = lambda wildcards: get_seed_msprime(wildcards.rep)
    shell:
        "python code/Simulate_Genotypes/generate_genotypes_SimpleGrid.py \
	       --outpre output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
	       --chr {params.chr_num} \
	       --sample_size 160 \
	       --length 10000 \
	       --Ne 1000 \
	       --mu 1e-07 \
	       --rho 1e-07 \
	       --tmove -9 \
	       --migrate 0.01 \
	       --seed {params.seed}"

rule format_VCF:
    input:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos_{chr}.vcf"
    output:
        gz="output/Simulate_Genotypes/SimpleGrid/{rep}/genos_{chr}.ids.vcf.gz"
    shell:
        """
	      head -n6 {input} > output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/header_{wildcards.chr}.txt
	            cat output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/header_{wildcards.chr}.txt <(cat output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos_{wildcards.chr}.vcf | awk -v OFS="\t" 'NR>6 {{$3=$1"_"$2"_A_T";$4="A"; $5="T"; print ;}}') | bgzip > {output.gz}
			        rm output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/header_{wildcards.chr}.txt
				      """

rule concat_vcfs:
    input:
        expand("output/Simulate_Genotypes/SimpleGrid/{{rep}}/genos_{chr}.ids.vcf.gz", chr=CHR)
    output:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.ids.vcf.gz"
    shell:
        """
        bcftools concat {input} -o output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/temp.vcf.gz -O z
        bcftools annotate --rename-chrs code/Simulate_Genotypes/convert_chr.txt output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/temp.vcf.gz -o {output} -O z
        rm output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/temp.vcf.gz
        """

rule convert_vcf_to_plink:
    input:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.ids.vcf.gz"
    output:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pgen",
      	"output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pvar"
    shell:
        "plink2 \
        --double-id \
        --make-pgen \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
        --vcf {input}"

rule create_panels_SimpleGrid:
    input:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pop"
    output:
        gwas_north="output/Simulate_Genotypes/SimpleGrid/{rep}/North/ids.gwas",
        gwas_south="output/Simulate_Genotypes/SimpleGrid/{rep}/South/ids.gwas",
        joint_gwas= "output/Simulate_Genotypes/SimpleGrid/{rep}/ids.joint_gwas",
	      test="output/Simulate_Genotypes/SimpleGrid/{rep}/ids.test"
    params:
        ss_test = SS_TEST
    shell:
        "Rscript code/Simulate_Genotypes/split_north-south_SimpleGrid.R {params.ss_test} {input} {output.gwas_north} {output.gwas_south} {output.joint_gwas} {output.test}"

rule split_into_test_gwas:
    input:
        gwas_north="output/Simulate_Genotypes/SimpleGrid/{rep}/North/ids.gwas",
	      test="output/Simulate_Genotypes/SimpleGrid/{rep}/ids.test",
	      gwas_south="output/Simulate_Genotypes/SimpleGrid/{rep}/South/ids.gwas",
	      gwas_joint="output/Simulate_Genotypes/SimpleGrid/{rep}/ids.joint_gwas",
	      psam="output/Simulate_Genotypes/SimpleGrid/{rep}/genos.psam",
		    pvar="output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pvar",
			  pgen="output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pgen"
    output:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.pgen",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.pvar",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.pgen",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.pvar",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.pgen",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.pvar",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.pgen",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.pvar"
    shell:
        """
	      plink2 \
	      --keep {input.gwas_north} \
		    --make-pgen \
			  --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas \
				--pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
				--rm-dup exclude-all

				plink2 \
	      --keep {input.gwas_south} \
		    --make-pgen \
			  --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas \
				--pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
				--rm-dup exclude-all

        plink2 \
        --keep {input.test} \
        --make-pgen \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
	      --rm-dup exclude-all

	      plink2 \
        --keep {input.gwas_joint} \
        --make-pgen \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos \
	      --rm-dup exclude-all
	      """

rule get_variant_freq:
    input:
          "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.psam",
	        "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.psam",
		      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.psam"
    output:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.afreq"
    shell:
        """
        plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas \
	      --freq \
		    --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas

			  plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas \
        --freq \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas

        plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test \
        --freq \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test
        """

rule get_common_snp_list:
    input:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.afreq",
    output:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/common_snp_ids.txt"
    script:
        "code/Simulate_Genotypes/get_common_snp_list_meta.R"

rule remake_panels_with_common_snps:
    input:
        common_id="output/Simulate_Genotypes/SimpleGrid/{rep}/common_snp_ids.txt",
	      test="output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.psam",
	      gwas_north="output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.psam",
        gwas_south="output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.psam",
        gwas_joint="output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.psam"
    output:
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test_common.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas_common.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas_common.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.psam",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test_common.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas_common.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas_common.afreq",
	      "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.afreq"
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas_common \
	      --make-pgen

        plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas_common \
	      --make-pgen

        plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test_common \
	      --make-pgen

	      plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas_common \
	      --make-pgen

	      plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas_common \
	      --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas_common \
		    --freq

		    plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas_common \
	      --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas_common \
		    --freq

        plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test_common \
	      --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-test_common \
		    --freq

		    plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas_common \
	      --out output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas_common \
		    --freq
        """

rule aggregate_genotypes:
    input:
        genos=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos_{chr}.vcf", chr=CHR, rep=REP),
        gz_chr=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos_{chr}.ids.vcf.gz", chr=CHR, rep=REP),
        frq_test=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.afreq", rep=REP),
        frq_gwas_north=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.afreq", rep=REP),
        frq_gwas_south=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.afreq", rep=REP),
        gz=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos.ids.vcf.gz", rep=REP),
	      gwas_pgen_north=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.pgen", rep=REP),
	      gwas_pvar_north=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.pvar", rep=REP),
	      gwas_psam_north=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas.psam", rep=REP),
	      gwas_pgen_joint=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.pgen", rep=REP),
	      gwas_pvar_joint=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.pvar", rep=REP),
	      gwas_psam_joint=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas.psam", rep=REP),
	      gwas_pgen_south=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.pgen", rep=REP),
	      gwas_pvar_south=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.pvar", rep=REP),
	      gwas_psam_south=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas.psam", rep=REP),
	      test_pgen=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.pgen", rep=REP),
	      test_pvar=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.pvar", rep=REP),
	      test_psam=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test.psam", rep=REP),
	      pgen=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pgen", rep=REP),
	      pvar=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pvar", rep=REP),
	      psam=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos.psam", rep=REP),
	      keep=expand("output/Simulate_Genotypes/SimpleGrid/{rep}/genos-test_common.psam", rep=REP)
    output:
        expand("output/Simulate_Genotypes/SimpleGrid/{rep}/ff.txt", rep=REP)
    shell:
        """
	 touch {output}
	 touch {input.keep}
	 rm {input.genos}
	 rm {input.frq_test}
	 rm {input.frq_gwas_north}
	 rm {input.frq_gwas_south}
	 rm {input.gz}
	 rm {input.gz_chr}
	 rm {input.gwas_pgen_north}
	 rm {input.gwas_pvar_north}
	 rm {input.gwas_psam_north}
	 rm {input.gwas_pgen_joint}
	 rm {input.gwas_pvar_joint}
	 rm {input.gwas_psam_joint}
	 rm {input.gwas_pgen_south}
	 rm {input.gwas_pvar_south}
	 rm {input.gwas_psam_south}
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
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.afreq"
    output:
        "output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.effects.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        seed = lambda wildcards: get_seed(wildcards.rep, wildcards.h2, wildcards.pheno, wildcards.env)
    shell:
        "Rscript code/Simulate_Phenotypes/draw_effect_sizes.R {input} {output} {params.her} 0.4 {params.seed}"

rule generate_genetic_values:
    input:
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.psam",
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.pvar",
        "output/Simulate_Genotypes/SimpleGrid/{rep}/genos-joint_gwas_common.pgen",
        "output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.effects.txt"
    output:
        "output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.gvalue.sscore"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/genos-joint_gwas_common \
	      --out output/Simulate_Phenotypes/SimpleGrid/{wildcards.rep}/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common.gvalue \
		    --score output/Simulate_Phenotypes/SimpleGrid/{wildcards.rep}/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common.effects.txt cols=dosagesum,scoresums"

rule simulate_phenotype_SimpleGrid:
    input:
        gvalues="output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.gvalue.sscore",
        pops="output/Simulate_Genotypes/SimpleGrid/{rep}/genos.pop"
    output:
        "output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.phenos.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        en = lambda wildcards: get_params(wildcards.env),
        seed = lambda wildcards: get_seed(wildcards.rep,wildcards.h2, wildcards.pheno, wildcards.env)
    shell:
        "Rscript code/Simulate_Phenotypes/simulate_phenotypes_SimpleGrid.R {input.gvalues} {input.pops} {output} {params.her} {params.en} {params.seed} {wildcards.pheno}"

# Run GWAS

rule gwas_North:
    input:
        genos="output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/SimpleGrid/{rep}/North/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.phenos.txt"
    output:
        "output/Run_GWAS/SimpleGrid/{rep}/North/{h2}/{pheno}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/North/genos-gwas_common \
        --read-freq {input.freq} \
        --glm allow-no-covars \
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/SimpleGrid/{wildcards.rep}/North/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common"

rule gwas_South:
    input:
        genos="output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/SimpleGrid/{rep}/South/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/SimpleGrid/{rep}/{h2}/{pheno}/{env}/genos-gwas_common.phenos.txt"
    output:
        "output/Run_GWAS/SimpleGrid/{rep}/South/{h2}/{pheno}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/South/genos-gwas_common \
        --read-freq {input.freq} \
        --glm allow-no-covars \
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/SimpleGrid/{wildcards.rep}/South/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common"

# FIX Ascertain SNPs for PRS

rule pick_SNPS:
    input:
        causal_effect="output/Simulate_Phenotypes/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    output:
        "output/PRS/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.c.betas",
        "output/PRS/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.c.p.betas",
        "output/PRS/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump.R {input.causal_effect} output/Run_GWAS/SimpleGrid/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common {params.pt} output/PRS/SimpleGrid/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-gwas_common"


# FIX Calculate true genetic value in test panel

rule calc_true_gv:
    input:
        genos="output/Simulate_Genotypes/SimpleGrid/{rep}/{config}/genos-test_common.psam",
        causal_effect="output/Simulate_Phenotypes/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-gwas_common.effects.txt",
        freq="output/Simulate_Genotypes/SimpleGrid/{rep}/{config}/genos-test_common.afreq"
    output:
        "output/PRS/SimpleGrid/{rep}/{config}/{h2}/{pheno}/{env}/genos-test_common.true.sscore",
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/SimpleGrid/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.causal_effect} cols=dosagesum,scoresums \
        --out output/PRS/SimpleGrid/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.pheno}/{wildcards.env}/genos-test_common.true \
        """




