#!/bin/sh
# properties = {"type": "single", "rule": "calc_prs", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.psam", "output/PRS/4PopSplit/V7/C2/genos-gwas_common.c.betas", "output/PRS/4PopSplit/V7/C2/genos-gwas_common.c.p.betas", "output/PRS/4PopSplit/V7/C2/genos-gwas_common.nc.betas", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.afreq"], "output": ["output/PRS/4PopSplit/V7/C2/genos-test_common.c.sscore", "output/PRS/4PopSplit/V7/C2/genos-test_common.c.p.sscore", "output/PRS/4PopSplit/V7/C2/genos-test_common.nc.sscore"], "wildcards": {"model": "4PopSplit", "rep": "V7", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 12, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/calc_prs.config=C2,model=4PopSplit,rep=V7.out", "error": "logs/calc_prs.config=C2,model=4PopSplit,rep=V7.err", "job-name": "calc_prs.config=C2,model=4PopSplit,rep=V7"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/PRS/4PopSplit/V7/C2/genos-test_common.c.p.sscore --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.psam output/PRS/4PopSplit/V7/C2/genos-gwas_common.c.betas output/PRS/4PopSplit/V7/C2/genos-gwas_common.c.p.betas output/PRS/4PopSplit/V7/C2/genos-gwas_common.nc.betas output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.afreq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules calc_prs --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/12.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/12.jobfailed; exit 1)

