#!/bin/sh
# properties = {"type": "single", "rule": "simulate_genotypes_4popsplit", "local": false, "input": [], "output": ["output/Simulate_Genotypes/4PopSplit/V1/genos_0.vcf", "output/Simulate_Genotypes/4PopSplit/V1/genos_1.vcf", "output/Simulate_Genotypes/4PopSplit/V1/genos.pop"], "wildcards": {"rep": "V1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 27, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/simulate_genotypes_4popsplit.rep=V1.out", "error": "logs/simulate_genotypes_4popsplit.rep=V1.err", "job-name": "simulate_genotypes_4popsplit.rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V1/genos.pop --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.prgff7m4 --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules simulate_genotypes_4popsplit --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.prgff7m4/27.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.prgff7m4/27.jobfailed; exit 1)

