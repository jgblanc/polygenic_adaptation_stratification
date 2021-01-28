#!/bin/sh
# properties = {"type": "single", "rule": "test_snp_freq", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.psam", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.pvar", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.pgen"], "output": ["output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.afreq"], "wildcards": {"model": "4PopSplit", "rep": "V1", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 10, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/test_snp_freq.config=C2,model=4PopSplit,rep=V1.out", "error": "logs/test_snp_freq.config=C2,model=4PopSplit,rep=V1.err", "job-name": "test_snp_freq.config=C2,model=4PopSplit,rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.afreq --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2bnf0mpv output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.psam output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.pvar output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.pgen --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules test_snp_freq --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2bnf0mpv/10.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.2bnf0mpv/10.jobfailed; exit 1)

