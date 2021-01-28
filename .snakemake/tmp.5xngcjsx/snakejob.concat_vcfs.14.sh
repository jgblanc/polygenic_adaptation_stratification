#!/bin/sh
# properties = {"type": "single", "rule": "concat_vcfs", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/genos_0.ids.vcf.gz", "output/Simulate_Genotypes/4PopSplit/V1/genos_1.ids.vcf.gz"], "output": ["output/Simulate_Genotypes/4PopSplit/V1/genos.ids.vcf.gz"], "wildcards": {"model": "4PopSplit", "rep": "V1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 14, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/concat_vcfs.model=4PopSplit,rep=V1.out", "error": "logs/concat_vcfs.model=4PopSplit,rep=V1.err", "job-name": "concat_vcfs.model=4PopSplit,rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V1/genos.ids.vcf.gz --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5xngcjsx output/Simulate_Genotypes/4PopSplit/V1/genos_0.ids.vcf.gz output/Simulate_Genotypes/4PopSplit/V1/genos_1.ids.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules concat_vcfs --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5xngcjsx/14.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5xngcjsx/14.jobfailed; exit 1)

