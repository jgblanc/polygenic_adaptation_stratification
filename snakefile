CHR =[]
for i in range(0, 200):
  CHR.append(str(i))
CONFIG=["C1"]
REP = []
for i in range(1,101):
  REP.append("A"+str(i))
HERITABILITY = ["h2-0.0"]
ENV=["env_0.0","env_0.5"]
#TS=["p-0.54", "p-0.57", "p-0.60","p-0.63","p-0.66"]
TS=["p-0.50"]
SIZE=2000
NUM_RESAMPLE=1000
PVALUE_THRESHOLD=1

wildcard_constraints:
    rep="[A-Z]\d+",
    config="C.",
    h2="h2-[0-1].[0-9]",
    env="env_-?[0-9].[0-9]*",
    ts="p-[0-1].[0-9][0-9]",
    dir="[a-z]*"

def get_seed_msprime(rep):
  out = int(''.join(list(rep)[1::])) * 1000
  return out

def get_params(x):
  out = x.split("-")[1]
  return out

def get_env(x):
  out = x.split("_")[1]
  return out

def get_seed_es(rep,config, h2, ts):
  rep = list(rep)[1]
  config = list(config)[1]
  h2 = h2.split(".")[1]
  ts_list = ts.split("-")[1].split(".")[1]
  ts_list = [int(i) for i in ts_list]
  ts = sum(ts_list)
  out = rep + config + h2  + str(ts)
  return out

def get_seed(rep,config, h2, ts, env):
  rep = list(rep)[1]
  config = list(config)[1]
  h2 = h2.split(".")[1]
  env = float(env.split("_")[1]) + 3.1415
  ts_list = ts.split("-")[1].split(".")[1]
  ts_list = [int(i) for i in ts_list]
  ts = sum(ts_list)
  out = int(env * float(rep + config + h2 + str(ts)))
  out = str(out)
  return out


rule all:
    input:
        expand("output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Qx.txt",rep=REP, config=CONFIG, h2=HERITABILITY, ts=TS, env=ENV)

# Simluate Genotypes

rule simulate_genotypes_4popsplit:
    output:
        expand("output/Simulate_Genotypes/4PopSplit/{{rep}}/genos_{chr}.vcf", chr=CHR),
	      "output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    params:
      chr_num = len(CHR),
      seed = lambda wildcards: get_seed_msprime(wildcards.rep)
    shell:
        "python code/Simulate_Genotypes/generate_genotypes_4PopSplit.py \
	       --outpre output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
	       --chr {params.chr_num} \
       	       --Nanc 10000 \
	       --NA 10000 \
	       --NB 10000 \
	       --NC 10000 \
	       --ND 10000 \
  	     -a 200 \
	       -b 200 \
	       -c 200 \
	       -d 200 \
         -s1 4400 \
          -s2 2200 \
          -L 10000"

rule format_VCF:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/genos_{chr}.vcf"
    output:
        gz="output/Simulate_Genotypes/4PopSplit/{rep}/genos_{chr}.ids.vcf.gz"
	      #csi="output/Simulate_Genotypes/4PopSplit/{rep}/genos_{chr}.ids.vcf.gz.csi"
    shell:
        """
	      head -n6 {input} > output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/header_{wildcards.chr}.txt
	            cat output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/header_{wildcards.chr}.txt <(cat output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos_{wildcards.chr}.vcf | awk -v OFS="\t" 'NR>6 {{$3=$1"_"$2"_A_T";$4="A"; $5="T"; print ;}}') | bgzip > {output.gz}
		          #bcftools index {output.gz}
			        rm output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/header_{wildcards.chr}.txt
				      """

rule concat_vcfs:
    input:
        expand("output/Simulate_Genotypes/4PopSplit/{{rep}}/genos_{chr}.ids.vcf.gz", chr=CHR)
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/genos.ids.vcf.gz"
    shell:
        """
        bcftools concat {input} -o output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/temp.vcf.gz -O z
        bcftools annotate --rename-chrs code/Simulate_Genotypes/convert_chr.txt output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/temp.vcf.gz -o {output} -O z
        rm output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/temp.vcf.gz
        """

rule convert_vcf_to_plink:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/genos.ids.vcf.gz"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/genos.psam",
	"output/Simulate_Genotypes/4PopSplit/{rep}/genos.pgen",
      	"output/Simulate_Genotypes/4PopSplit/{rep}/genos.pvar"
    shell:
        "plink2 \
        --double-id \
        --make-pgen \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
        --vcf {input}"

rule create_panels_4PopSplit:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/C1/ids.gwas",
	      "output/Simulate_Genotypes/4PopSplit/{rep}/C1/ids.test",
	            "output/Simulate_Genotypes/4PopSplit/{rep}/C2/ids.gwas",
		          "output/Simulate_Genotypes/4PopSplit/{rep}/C2/ids.test"
    script:
        "code/Simulate_Genotypes/split_gwas-test_4PopSplit.R"

rule split_into_test_gwas:
    input:
        gwas="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/ids.gwas",
	      test="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/ids.test",
	            psam="output/Simulate_Genotypes/4PopSplit/{rep}/genos.psam",
		          pvar="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pvar",
			        pgen="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pgen"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.psam",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pgen",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pvar",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.psam",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pgen",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pvar"
    shell:
        """
	      plink2 \
	            --keep {input.gwas} \
		          --make-pgen \
			        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas \
				      --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
				            --rm-dup exclude-all

        plink2 \
        --keep {input.test} \
        --make-pgen \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test-big \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/genos \
	      --rm-dup exclude-all
	            """

rule downsample_test:
    input:
      psam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.psam",
          pgen="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pgen",
	      pvar="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pvar"
    output:
      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.psam",
      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pvar",
      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pgen",
      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/downsample.id"
    params:
      size = SIZE
    shell:
        """
	set +o pipefail;
        awk 'NR > 1' {input.psam} | cut -f 1,2 | sort -R | head -n {params.size} > output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/downsample.id

        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test-big \
        --keep output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/downsample.id \
        --make-pgen \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test
        """

rule get_variant_freq:
    input:
          "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.psam",
	        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pvar",
		      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pgen"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.afreq",
	"output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.afreq"
    shell:
        """
        plink2 \
	      --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test \
	            --freq \
		          --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test

			        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas \
        --freq \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas
        """

rule get_common_snp_list:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.afreq",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.afreq"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/common_snp_ids.txt"
    script:
        "code/Simulate_Genotypes/get_common_snp_list.R"

rule remake_panels_with_common_snps:
    input:
        common_id="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/common_snp_ids.txt",
	      test_psam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.psam",
	            test_pvar="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pvar",
		          test_pgen="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pgen",
			        gwas_psam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.psam",
				      gwas_pvar="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pvar",
				            gwas_pgen="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pgen"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pvar",
	      "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pgen",
	            "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
		          "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pvar",
			        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pgen"
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
	      --make-pgen

        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas \
        --extract {input.common_id} \
        --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	      --make-pgen
        """

rule common_snp_freq:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pgen"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	            --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
		          --freq"

rule aggregate_genotypes:
    input:
        frq=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq", rep=REP, config=CONFIG),
        genos=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos_{chr}.vcf", chr=CHR, rep=REP, config=CONFIG),
        gz_chr=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos_{chr}.ids.vcf.gz", chr=CHR, rep=REP),
        frq_test=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.afreq", rep=REP, config=CONFIG),
        frq_gwas=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.afreq", rep=REP, config=CONFIG),
        gz=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos.ids.vcf.gz", rep=REP),
	      gwas_pgen=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pgen", rep=REP, config=CONFIG),
	      gwas_pvar=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.pvar", rep=REP, config=CONFIG),
	      gwas_psam=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas.psam", rep=REP, config=CONFIG),
	      test_pgen=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pgen", rep=REP, config=CONFIG),
	      test_pvar=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.pvar", rep=REP, config=CONFIG),
	      test_psam=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test.psam", rep=REP, config=CONFIG),
	      pgen=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos.pgen", rep=REP, config=CONFIG),
	      pvar=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos.pvar", rep=REP, config=CONFIG),
	      psam=expand("output/Simulate_Genotypes/4PopSplit/{rep}/genos.psam", rep=REP, config=CONFIG),
	      big_psam=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.psam", rep=REP, config=CONFIG),
	      big_pvar=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pvar", rep=REP, config=CONFIG),
	      big_pgen=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test-big.pgen", rep=REP, config=CONFIG),
        id=expand("output/Simulate_Genotypes/4PopSplit/{rep}/{config}/downsample.id", rep=REP, config=CONFIG)
    output:
        expand("output/Simulate_Genotypes/4PopSplit/{rep}/ff.txt", rep=REP)
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
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        "output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        seed = lambda wildcards: get_seed_es(wildcards.rep, wildcards.config, wildcards.h2, wildcards.ts),
        prob = lambda wildcards: get_params(wildcards.ts)
    shell:
        "Rscript code/Simulate_Phenotypes/draw_effect_sizes_4PopSplit.R {input.freq} {output} {params.her} 0.4 {params.seed} output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common {input.pops} {params.prob}"


rule generate_genetic_values:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pvar",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.pgen",
        "output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt"
    output:
        "output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.gvalue.sscore"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
	            --out output/Simulate_Phenotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/genos-gwas_common.gvalue \
		          --score output/Simulate_Phenotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/genos-gwas_common.effects.txt cols=dosagesum,scoresums"


rule simulate_phenotype_4PopSplit:
    input:
        gvalues="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.gvalue.sscore",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop"
    output:
        "output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.phenos.txt"
    params:
        her = lambda wildcards: get_params(wildcards.h2),
        en = lambda wildcards: get_env(wildcards.env),
        seed = lambda wildcards: get_seed(wildcards.rep, wildcards.config, wildcards.h2, wildcards.ts, wildcards.env)
    shell:
        "Rscript code/Simulate_Phenotypes/simulate_phenotypes_4PopSplit.R {input.gvalues} {input.pops} {output} {params.her} {params.en} {params.seed}"


# Run GWAS

rule gwas_no_correction:
    input:
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.phenos.txt"
    output:
        "output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm allow-no-covars\
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common"



# Ascertain SNPs to include in PRS

rule pick_SNPS:
    input:
        causal_effect="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.pheno_strat.glm.linear"
    output:
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.p.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump.R {input.causal_effect} output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common {params.pt} output/PRS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common"

# Get allele freq in test panel

rule test_snp_freq:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pvar",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pgen"
    output:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq"
    shell:
        "plink2 \
	      --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
	            --out output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
		          --freq"

# Calculate true breeding value

rule calc_true_gv:
    input:
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        causal_effect="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt",
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq"
    output:
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-test_common.true.sscore",
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.causal_effect} cols=dosagesum,scoresums \
        --out output/PRS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/genos-test_common.true \
        """

## Include TGWAS as a covariate

# Generate Test Vector

rule make_test_vector:
    input:
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        fam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam"
    output:
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/Tvec.txt"
    shell:
        "Rscript code/Calculate_TGWAS/4PopSplit_make_tvec.R {input.pops} {input.fam} {output}"

# Project T using CGD

rule calc_TGWAS:
    input:
        genos_test="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        genos_gwas="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        Tvec="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/Tvec.txt"
    output:
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/TGWAS.txt"
    shell:
        """
        Rscript code/Calculate_TGWAS/Compute_TGWAS.R output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common {input.Tvec} output/Calculate_TGWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/
        rm output/Calculate_TGWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/[Gxpb]*
	      """

# Format Covariate file

rule format_covars:
    input:
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        fam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        TGWAS="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/TGWAS.txt"
    output:
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/TGWAS-ID_covars.txt"
    shell:
        """
        Rscript code/Calculate_TGWAS/format_ID_covars.R {input.pops} {input.TGWAS} {input.fam} {output}
        """


# Re-run GWAS

rule gwas_TGWAS:
    input:
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.phenos.txt",
        TGWAS="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/TGWAS-ID_covars.txt"
    output:
        "output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm hide-covar \
        --covar {input.TGWAS} \
        --covar-col-nums 3 \
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-TGWAS"

rule gwas_PopID:
    input:
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.psam",
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-gwas_common.afreq",
        pheno="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.phenos.txt",
        TGWAS="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/TGWAS-ID_covars.txt"
    output:
        "output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.pheno_strat.glm.linear"
    shell:
        "plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-gwas_common \
        --read-freq {input.freq} \
        --glm hide-covar \
        --covar {input.TGWAS} \
        --covar-col-nums 4 \
        --pheno {input.pheno} \
        --pheno-name pheno_strat \
        --out output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-ID"

# Remake PRS

rule pick_SNPS_TGWAS:
    input:
        causal_effect="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.pheno_strat.glm.linear"
    output:
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.c.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.c.p.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump.R {input.causal_effect} output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-TGWAS {params.pt} output/PRS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-TGWAS"

rule pick_SNPS_ID:
    input:
        causal_effect="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt",
        gwas_strat="output/Run_GWAS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.pheno_strat.glm.linear"
    output:
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.p.betas",
        "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.nc.betas"
    params:
        pt = PVALUE_THRESHOLD
    shell:
        "Rscript code/PRS/clump.R {input.causal_effect} output/Run_GWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-ID {params.pt} output/PRS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-gwas_common-ID"


# Do polygenic adaptation test

rule PCA_test:
    input:
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pgen",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.pvar",
        "output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/Tvec.txt"
    params:
        n_minus_1 = int(SIZE)-1,
        col_start = 6,
        col_end = int(SIZE) + 4
    output:
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/pca.eigenvec",
        "output/Calculate_TGWAS/4PopSplit/{rep}/{config}/pca.eigenval",
    shell:
        """
        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
       --pca {params.n_minus_1} \
       --out output/Calculate_TGWAS/4PopSplit/{wildcards.rep}/{wildcards.config}/pca
        """


rule calc_lambdaT:
    input:
        vecs="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/pca.eigenvec",
        vals="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/pca.eigenval",
        tvec="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/Tvec.txt",
    output:
        "output/PGA_test/4PopSplit/{rep}/{config}/Lambda_T.txt"
    shell:
        """
	Rscript code/PGA_test/4PopSplit_calc_lambdaT.R {input.vecs} {input.vals} {input.tvec} {output}
	"""

rule calc_Va:
    input:
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq",
        c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.betas",
        cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.p.betas",
        nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.nc.betas"
    output:
        "output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Va.txt"
    shell:
      """
          Rscript code/PGA_test/calc_Va_strat.R {input.freq} {input.c} {input.cp} {input.nc} {output}
	      """

rule calc_Va_Tm:
    input:
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq",
        c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.c.betas",
        cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.c.p.betas",
        nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.nc.betas"
    output:
        "output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Va-Tm.txt"
    shell:
      """
          Rscript code/PGA_test/calc_Va_strat.R {input.freq} {input.c} {input.cp} {input.nc} {output}
	      """

rule calc_Va_ID:
    input:
        freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq",
        c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.betas",
        cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.p.betas",
        nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.nc.betas"
    output:
        "output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Va-ID.txt"
    shell:
      """
          Rscript code/PGA_test/calc_Va_strat.R {input.freq} {input.c} {input.cp} {input.nc} {output}
	      """


rule Calc_Qx:
    input:
        c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.betas",
        cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.p.betas",
        nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.nc.betas",
        c_TGWAS="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.c.betas",
        cp_TGWAS="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.c.p.betas",
        nc_TGWAS="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-TGWAS.nc.betas",
        c_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.betas",
        cp_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.p.betas",
        nc_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.nc.betas",
        genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        lambda_T="output/PGA_test/4PopSplit/{rep}/{config}/Lambda_T.txt",
        true="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-test_common.true.sscore",
        Tvec="output/Calculate_TGWAS/4PopSplit/{rep}/{config}/Tvec.txt",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        es="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt"
    output:
        qx="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Qx.txt",
        pgs="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/PGS.txt"
    params:
        num=NUM_RESAMPLE
    shell:
      """
          Rscript code/PGA_test/calc_Qx_4PopSplit.R {input.c} {input.cp} {input.nc} {input.c_TGWAS} {input.cp_TGWAS} {input.nc_TGWAS} {input.c_ID} {input.cp_ID} {input.nc_ID} output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common {input.lambda_T} {input.true} {input.Tvec} {input.pops} {params.num} {output.qx} {output.pgs} {input.es}
	      """


# Calculate true signal magnitude

rule calc_ts_magnitude:
    input:
        psam="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
        pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
        true_es="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/genos-gwas_common.effects.txt"
    output:
        "output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/ts_magnitude.txt"
    shell:
        "Rscript code/Simulate_Phenotypes/calculate_true_signal_magnitude.R {input.pops} {input.true_es} output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common {output}"



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

