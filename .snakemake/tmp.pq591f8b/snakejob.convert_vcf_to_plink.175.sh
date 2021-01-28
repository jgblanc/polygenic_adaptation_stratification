#!/bin/sh
# properties = {"type": "single", "rule": "convert_vcf_to_plink", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V4/genos.ids.vcf.gz"], "output": ["output/Simulate_Genotypes/4PopSplit/V4/genos.psam", "output/Simulate_Genotypes/4PopSplit/V4/genos.pgen", "output/Simulate_Genotypes/4PopSplit/V4/genos.pvar"], "wildcards": {"model": "4PopSplit", "rep": "V4"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 175, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/convert_vcf_to_plink.model=4PopSplit,rep=V4.out", "error": "logs/convert_vcf_to_plink.model=4PopSplit,rep=V4.err", "job-name": "convert_vcf_to_plink.model=4PopSplit,rep=V4"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V4/genos.psam --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b output/Simulate_Genotypes/4PopSplit/V4/genos.ids.vcf.gz --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules convert_vcf_to_plink --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/175.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.pq591f8b/175.jobfailed; exit 1)

