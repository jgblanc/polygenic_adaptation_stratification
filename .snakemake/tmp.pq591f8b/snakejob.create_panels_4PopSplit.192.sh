#!/bin/sh
# properties = {"type": "single", "rule": "create_panels_4PopSplit", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V7/genos.pop"], "output": ["output/Simulate_Genotypes/4PopSplit/V7/C1/ids.gwas", "output/Simulate_Genotypes/4PopSplit/V7/C1/ids.test", "output/Simulate_Genotypes/4PopSplit/V7/C2/ids.gwas", "output/Simulate_Genotypes/4PopSplit/V7/C2/ids.test"], "wildcards": {"model": "4PopSplit", "rep": "V7"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 192, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/create_panels_4PopSplit.model=4PopSplit,rep=V7.out", "error": "logs/create_panels_4PopSplit.model=4PopSplit,rep=V7.err", "job-name": "create_panels_4PopSplit.model=4PopSplit,rep=V7"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V7/C1/ids.gwas --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Genotypes/4PopSplit/V7/genos.pop --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules create_panels_4PopSplit --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/192.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/192.jobfailed; exit 1)
