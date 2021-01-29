#!/bin/sh
# properties = {"type": "single", "rule": "get_variant_freq", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.psam", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pvar", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pgen"], "output": ["output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.afreq", "output/Simulate_Genotypes/4PopSplit/V2/C1/genos-gwas.afreq"], "wildcards": {"model": "4PopSplit", "rep": "V2", "config": "C1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 21, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/get_variant_freq.config=C1,model=4PopSplit,rep=V2.out", "error": "logs/get_variant_freq.config=C1,model=4PopSplit,rep=V2.err", "job-name": "get_variant_freq.config=C1,model=4PopSplit,rep=V2"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.afreq --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.psam output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pvar output/Simulate_Genotypes/4PopSplit/V2/C1/genos-test.pgen --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules get_variant_freq --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/21.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/21.jobfailed; exit 1)

