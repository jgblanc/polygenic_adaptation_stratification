#!/bin/sh
# properties = {"type": "single", "rule": "concat_vcfs", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V4/genos_0.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_1.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_2.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_3.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_4.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_5.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_6.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_7.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_8.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_9.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_10.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_11.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_12.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_13.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_14.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_15.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_16.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_17.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_18.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V4/genos_19.ids.vcf.gz"], "output": ["output/Simulate_Genotypes/4PopSplit/V4/genos.ids.vcf.gz"], "wildcards": {"model": "4PopSplit", "rep": "V4"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 188, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/concat_vcfs.model=4PopSplit,rep=V4.out", "error": "logs/concat_vcfs.model=4PopSplit,rep=V4.err", "job-name": "concat_vcfs.model=4PopSplit,rep=V4"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V4/genos.ids.vcf.gz --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Genotypes/4PopSplit/V4/genos_0.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_1.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_2.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_3.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_4.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_5.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_6.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_7.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_8.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_9.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_10.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_11.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_12.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_13.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_14.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_15.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_16.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_17.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_18.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V4/genos_19.ids.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules concat_vcfs --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/188.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/188.jobfailed; exit 1)

