#!/bin/sh
# properties = {"type": "single", "rule": "get_common_snp_list", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/C1/genos-test.afreq", "output/Simulate_Genotypes/4PopSplit/V1/C1/genos-gwas.afreq"], "output": ["output/Simulate_Genotypes/4PopSplit/V1/C1/common_snp_ids.txt"], "wildcards": {"model": "4PopSplit", "rep": "V1", "config": "C1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 13, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/get_common_snp_list.config=C1,model=4PopSplit,rep=V1.out", "error": "logs/get_common_snp_list.config=C1,model=4PopSplit,rep=V1.err", "job-name": "get_common_snp_list.config=C1,model=4PopSplit,rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V1/C1/common_snp_ids.txt --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.tstchr36 output/Simulate_Genotypes/4PopSplit/V1/C1/genos-test.afreq output/Simulate_Genotypes/4PopSplit/V1/C1/genos-gwas.afreq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules get_common_snp_list --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.tstchr36/13.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.tstchr36/13.jobfailed; exit 1)

