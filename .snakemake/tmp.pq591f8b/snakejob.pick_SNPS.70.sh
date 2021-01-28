#!/bin/sh
# properties = {"type": "single", "rule": "pick_SNPS", "local": false, "input": ["output/Simulate_Phenotypes/4PopSplit/V8/C1/genos-gwas_common.effects.txt", "output/Run_GWAS/4PopSplit/V8/C1/genos-gwas_common.pheno_random.glm.linear", "output/Run_GWAS/4PopSplit/V8/C1/genos-gwas_common.pheno_strat.glm.linear"], "output": ["output/PRS/4PopSplit/V8/C1/genos-gwas_common.c.betas", "output/PRS/4PopSplit/V8/C1/genos-gwas_common.c.p.betas", "output/PRS/4PopSplit/V8/C1/genos-gwas_common.nc.betas"], "wildcards": {"model": "4PopSplit", "rep": "V8", "config": "C1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 70, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/pick_SNPS.config=C1,model=4PopSplit,rep=V8.out", "error": "logs/pick_SNPS.config=C1,model=4PopSplit,rep=V8.err", "job-name": "pick_SNPS.config=C1,model=4PopSplit,rep=V8"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/PRS/4PopSplit/V8/C1/genos-gwas_common.c.betas --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Phenotypes/4PopSplit/V8/C1/genos-gwas_common.effects.txt output/Run_GWAS/4PopSplit/V8/C1/genos-gwas_common.pheno_random.glm.linear output/Run_GWAS/4PopSplit/V8/C1/genos-gwas_common.pheno_strat.glm.linear --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules pick_SNPS --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/70.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/70.jobfailed; exit 1)

