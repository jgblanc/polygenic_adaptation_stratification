#!/bin/sh
# properties = {"type": "single", "rule": "generate_genetic_values", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.psam", "output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.pvar", "output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.pgen", "output/Simulate_Phenotypes/4PopSplit/V6/C2/genos-gwas_common.effects.txt"], "output": ["output/Simulate_Phenotypes/4PopSplit/V6/C2/genos-gwas_common.gvalue.sscore"], "wildcards": {"model": "4PopSplit", "rep": "V6", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 228, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/generate_genetic_values.config=C2,model=4PopSplit,rep=V6.out", "error": "logs/generate_genetic_values.config=C2,model=4PopSplit,rep=V6.err", "job-name": "generate_genetic_values.config=C2,model=4PopSplit,rep=V6"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Phenotypes/4PopSplit/V6/C2/genos-gwas_common.gvalue.sscore --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.psam output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.pvar output/Simulate_Genotypes/4PopSplit/V6/C2/genos-gwas_common.pgen output/Simulate_Phenotypes/4PopSplit/V6/C2/genos-gwas_common.effects.txt --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules generate_genetic_values --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/228.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/228.jobfailed; exit 1)
