

# Ascertain SNPs to include in PRS



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
causal_effect="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.effects.txt",
freq="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.afreq"
output:
  "output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-test_common.true.sscore",
shell:
  """
        plink2 \
        --pfile output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common \
        --read-freq {input.freq} \
        --score {input.causal_effect} cols=dosagesum,scoresums \
        --out output/PRS/4PopSplit/{wildcards.rep}/{wildcards.config}/{wildcards.h2}/{wildcards.ts}/{wildcards.env}/genos-test_common.true \
        """

## Include Tm as a covariate






# Remake PRS






rule calc_lambdaT:
  input:
  vecs="output/Calculate_Tm/4PopSplit/{rep}/{config}/pca.eigenvec",
vals="output/Calculate_Tm/4PopSplit/{rep}/{config}/pca.eigenval",
tvec="output/Calculate_Tm/4PopSplit/{rep}/{config}/Tvec.txt",
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

rule Calc_Qx:
  input:
  c="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.betas",
cp="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.c.p.betas",
nc="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.nc.betas",
c_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.c.betas",
cp_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.c.p.betas",
nc_Tm="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-Tm.nc.betas",
c_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.betas",
cp_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.c.p.betas",
nc_ID="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common-ID.nc.betas",
genos="output/Simulate_Genotypes/4PopSplit/{rep}/{config}/genos-test_common.psam",
lambda_T="output/PGA_test/4PopSplit/{rep}/{config}/Lambda_T.txt",
true="output/PRS/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-test_common.true.sscore",
Tvec="output/Calculate_Tm/4PopSplit/{rep}/{config}/Tvec.txt",
pops="output/Simulate_Genotypes/4PopSplit/{rep}/genos.pop",
es="output/Simulate_Phenotypes/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/genos-gwas_common.effects.txt"
output:
  qx="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/Qx.txt",
pgs="output/PGA_test/4PopSplit/{rep}/{config}/{h2}/{ts}/{env}/PGS.txt"
params:
  num=NUM_RESAMPLE
shell:
  """
          Rscript code/PGA_test/calc_Qx_4PopSplit.R {input.c} {input.cp} {input.nc} {input.c_Tm} {input.cp_Tm} {input.nc_Tm} {input.c_ID} {input.cp_ID} {input.nc_ID} output/Simulate_Genotypes/4PopSplit/{wildcards.rep}/{wildcards.config}/genos-test_common {input.lambda_T} {input.true} {input.Tvec} {input.pops} {params.num} {output.qx} {output.pgs} {input.es}
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






