#!/bin/sh
# properties = {"type": "single", "rule": "aggregate_genotypes", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/C1/genos-gwas_common.afreq", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.afreq"], "output": [], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 0, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/aggregate_genotypes..out", "error": "logs/aggregate_genotypes..err", "job-name": "aggregate_genotypes."}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake aggregate_genotypes --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2mdnbm81 output/Simulate_Genotypes/4PopSplit/V1/C1/genos-gwas_common.afreq output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.afreq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules aggregate_genotypes --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2mdnbm81/0.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2mdnbm81/0.jobfailed; exit 1)

