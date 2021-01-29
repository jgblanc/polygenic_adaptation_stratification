#!/bin/sh
# properties = {"type": "single", "rule": "simulate_genotypes_4popsplit", "local": false, "input": [], "output": ["output/Simulate_Genotypes/4PopSplit/V6/genos_0.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_1.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_2.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_3.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_4.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_5.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_6.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_7.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_8.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_9.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_10.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_11.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_12.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_13.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_14.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_15.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_16.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_17.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_18.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos_19.vcf", "output/Simulate_Genotypes/4PopSplit/V6/genos.pop"], "wildcards": {"rep": "V6"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 195, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/simulate_genotypes_4popsplit.rep=V6.out", "error": "logs/simulate_genotypes_4popsplit.rep=V6.err", "job-name": "simulate_genotypes_4popsplit.rep=V6"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V6/genos_0.vcf --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulate_genotypes_4popsplit --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/195.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/195.jobfailed; exit 1)

