#!/bin/sh
# properties = {"type": "single", "rule": "remake_panels_with_common_snps", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V7/C2/common_snp_ids.txt", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.psam", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.pvar", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.pgen", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.psam", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.pvar", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.pgen"], "output": ["output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.psam", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.pvar", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test_common.pgen", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas_common.psam", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas_common.pvar", "output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas_common.pgen"], "wildcards": {"model": "4PopSplit", "rep": "V7", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 56, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/remake_panels_with_common_snps.config=C2,model=4PopSplit,rep=V7.out", "error": "logs/remake_panels_with_common_snps.config=C2,model=4PopSplit,rep=V7.err", "job-name": "remake_panels_with_common_snps.config=C2,model=4PopSplit,rep=V7"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas_common.psam --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Genotypes/4PopSplit/V7/C2/common_snp_ids.txt output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.psam output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.pvar output/Simulate_Genotypes/4PopSplit/V7/C2/genos-test.pgen output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.psam output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.pvar output/Simulate_Genotypes/4PopSplit/V7/C2/genos-gwas.pgen --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules remake_panels_with_common_snps --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/56.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/56.jobfailed; exit 1)

