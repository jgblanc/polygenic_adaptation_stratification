#!/bin/sh
# properties = {"type": "single", "rule": "common_snp_freq", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.psam", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.pvar", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.pgen"], "output": ["output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.afreq"], "wildcards": {"model": "4PopSplit", "rep": "V1", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 20, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/common_snp_freq.config=C2,model=4PopSplit,rep=V1.out", "error": "logs/common_snp_freq.config=C2,model=4PopSplit,rep=V1.err", "job-name": "common_snp_freq.config=C2,model=4PopSplit,rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.afreq --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5hqx_r3u output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.psam output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.pvar output/Simulate_Genotypes/4PopSplit/V1/C2/genos-gwas_common.pgen --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules common_snp_freq --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5hqx_r3u/20.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.5hqx_r3u/20.jobfailed; exit 1)
