#!/bin/sh
# properties = {"type": "single", "rule": "draw_effect_sizes", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V7/C1/genos-gwas_common.afreq"], "output": ["output/Simulate_Phenotypes/4PopSplit/V7/C1/genos-gwas_common.effects.txt"], "wildcards": {"model": "4PopSplit", "rep": "V7", "config": "C1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 91, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/draw_effect_sizes.config=C1,model=4PopSplit,rep=V7.out", "error": "logs/draw_effect_sizes.config=C1,model=4PopSplit,rep=V7.err", "job-name": "draw_effect_sizes.config=C1,model=4PopSplit,rep=V7"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Phenotypes/4PopSplit/V7/C1/genos-gwas_common.effects.txt --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Genotypes/4PopSplit/V7/C1/genos-gwas_common.afreq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules draw_effect_sizes --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/91.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/91.jobfailed; exit 1)
