#!/bin/sh
# properties = {"type": "single", "rule": "aggregate_genotypes", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas_common.afreq", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas_common.afreq", "output/Simulate_Genotypes/4PopSplit/V2/genos_0.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_1.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_2.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_3.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_4.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_5.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_6.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_7.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_8.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_9.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_10.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_11.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_12.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_13.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_14.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_15.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_16.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_17.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_18.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_19.vcf", "output/Simulate_Genotypes/4PopSplit/V2/genos_0.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_1.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_2.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_3.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_4.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_5.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_6.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_7.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_8.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_9.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_10.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_11.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_12.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_13.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_14.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_15.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_16.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_17.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_18.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/genos_19.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.afreq", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.afreq", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.afreq", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.afreq", "output/Simulate_Genotypes/4PopSplit/V2/genos.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.pgen", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.pgen", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.pvar", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.pvar", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.psam", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.psam", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pgen", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.pgen", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pvar", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.pvar", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.psam", "output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.psam", "output/Simulate_Genotypes/4PopSplit/V2/genos.pgen", "output/Simulate_Genotypes/4PopSplit/V2/genos.pvar", "output/Simulate_Genotypes/4PopSplit/V2/genos.psam"], "output": ["output/Simulate_Genotypes/4PopSplit/V2/ff.txt"], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 51, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/aggregate_genotypes..out", "error": "logs/aggregate_genotypes..err", "job-name": "aggregate_genotypes."}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V2/ff.txt --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas_common.afreq output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas_common.afreq output/Simulate_Genotypes/4PopSplit/V2/genos_0.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_1.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_2.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_3.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_4.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_5.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_6.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_7.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_8.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_9.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_10.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_11.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_12.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_13.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_14.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_15.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_16.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_17.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_18.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_19.vcf output/Simulate_Genotypes/4PopSplit/V2/genos_0.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_1.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_2.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_3.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_4.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_5.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_6.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_7.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_8.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_9.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_10.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_11.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_12.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_13.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_14.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_15.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_16.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_17.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_18.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/genos_19.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.afreq output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.afreq output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.afreq output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.afreq output/Simulate_Genotypes/4PopSplit/V2/genos.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.pgen output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.pgen output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.pvar output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.pvar output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.psam output/Simulate_Genotypes/4PopSplit/V2/C2/genos-gwas.psam output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pgen output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.pgen output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pvar output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.pvar output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.psam output/Simulate_Genotypes/4PopSplit/V2/C2/genos-test.psam output/Simulate_Genotypes/4PopSplit/V2/genos.pgen output/Simulate_Genotypes/4PopSplit/V2/genos.pvar output/Simulate_Genotypes/4PopSplit/V2/genos.psam --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules aggregate_genotypes --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/51.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/51.jobfailed; exit 1)

