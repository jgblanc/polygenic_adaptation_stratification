#!/bin/sh
# properties = {"type": "single", "rule": "split_into_test_gwas", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V9/C1/ids.gwas", "output/Simulate_Genotypes/4PopSplit/V9/C1/ids.test", "output/Simulate_Genotypes/4PopSplit/V9/genos.psam", "output/Simulate_Genotypes/4PopSplit/V9/genos.pvar", "output/Simulate_Genotypes/4PopSplit/V9/genos.pgen"], "output": ["output/Simulate_Genotypes/4PopSplit/V9/C1/genos-test.psam", "output/Simulate_Genotypes/4PopSplit/V9/C1/genos-test.pgen", "output/Simulate_Genotypes/4PopSplit/V9/C1/genos-test.pvar", "output/Simulate_Genotypes/4PopSplit/V9/C1/genos-gwas.psam", "output/Simulate_Genotypes/4PopSplit/V9/C1/genos-gwas.pgen", "output/Simulate_Genotypes/4PopSplit/V9/C1/genos-gwas.pvar"], "wildcards": {"model": "4PopSplit", "rep": "V9", "config": "C1"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 122, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/split_into_test_gwas.config=C1,model=4PopSplit,rep=V9.out", "error": "logs/split_into_test_gwas.config=C1,model=4PopSplit,rep=V9.err", "job-name": "split_into_test_gwas.config=C1,model=4PopSplit,rep=V9"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V9/C1/genos-gwas.pgen --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Genotypes/4PopSplit/V9/C1/ids.gwas output/Simulate_Genotypes/4PopSplit/V9/C1/ids.test output/Simulate_Genotypes/4PopSplit/V9/genos.psam output/Simulate_Genotypes/4PopSplit/V9/genos.pvar output/Simulate_Genotypes/4PopSplit/V9/genos.pgen --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules split_into_test_gwas --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/122.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/122.jobfailed; exit 1)

