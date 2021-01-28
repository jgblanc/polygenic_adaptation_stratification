#!/bin/sh
# properties = {"type": "single", "rule": "calc_true_gv", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.psam", "output/Simulate_Phenotypes/4PopSplit/V1/C2/genos-gwas_common.effects.txt", "output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.afreq"], "output": ["output/PRS/4PopSplit/V1/C2/genos-test_common.true.sscore"], "wildcards": {"model": "4PopSplit", "rep": "V1", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 4, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/calc_true_gv.config=C2,model=4PopSplit,rep=V1.out", "error": "logs/calc_true_gv.config=C2,model=4PopSplit,rep=V1.err", "job-name": "calc_true_gv.config=C2,model=4PopSplit,rep=V1"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/PRS/4PopSplit/V1/C2/genos-test_common.true.sscore --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.m82vowb0 output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.psam output/Simulate_Phenotypes/4PopSplit/V1/C2/genos-gwas_common.effects.txt output/Simulate_Genotypes/4PopSplit/V1/C2/genos-test_common.afreq --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules calc_true_gv --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.m82vowb0/4.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.m82vowb0/4.jobfailed; exit 1)

