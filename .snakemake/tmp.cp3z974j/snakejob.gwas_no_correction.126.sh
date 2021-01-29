#!/bin/sh
# properties = {"type": "single", "rule": "gwas_no_correction", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V9/C2/genos-gwas_common.psam", "output/Simulate_Genotypes/4PopSplit/V9/C2/genos-gwas_common.afreq", "output/Simulate_Phenotypes/4PopSplit/V9/C2/genos-gwas_common.phenos.txt"], "output": ["output/Run_GWAS/4PopSplit/V9/C2/genos-gwas_common.pheno_random.glm.linear", "output/Run_GWAS/4PopSplit/V9/C2/genos-gwas_common.pheno_strat.glm.linear"], "wildcards": {"model": "4PopSplit", "rep": "V9", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 126, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/gwas_no_correction.config=C2,model=4PopSplit,rep=V9.out", "error": "logs/gwas_no_correction.config=C2,model=4PopSplit,rep=V9.err", "job-name": "gwas_no_correction.config=C2,model=4PopSplit,rep=V9"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Run_GWAS/4PopSplit/V9/C2/genos-gwas_common.pheno_random.glm.linear --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Genotypes/4PopSplit/V9/C2/genos-gwas_common.psam output/Simulate_Genotypes/4PopSplit/V9/C2/genos-gwas_common.afreq output/Simulate_Phenotypes/4PopSplit/V9/C2/genos-gwas_common.phenos.txt --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules gwas_no_correction --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/126.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/126.jobfailed; exit 1)

