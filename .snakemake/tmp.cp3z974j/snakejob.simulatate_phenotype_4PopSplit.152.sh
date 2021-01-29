#!/bin/sh
# properties = {"type": "single", "rule": "simulatate_phenotype_4PopSplit", "local": false, "input": ["output/Simulate_Phenotypes/4PopSplit/V4/C2/genos-gwas_common.gvalue.sscore", "output/Simulate_Genotypes/4PopSplit/V4/genos.pop"], "output": ["output/Simulate_Phenotypes/4PopSplit/V4/C2/genos-gwas_common.phenos.txt"], "wildcards": {"model": "4PopSplit", "rep": "V4", "config": "C2"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 152, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/simulatate_phenotype_4PopSplit.config=C2,model=4PopSplit,rep=V4.out", "error": "logs/simulatate_phenotype_4PopSplit.config=C2,model=4PopSplit,rep=V4.err", "job-name": "simulatate_phenotype_4PopSplit.config=C2,model=4PopSplit,rep=V4"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Phenotypes/4PopSplit/V4/C2/genos-gwas_common.phenos.txt --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Phenotypes/4PopSplit/V4/C2/genos-gwas_common.gvalue.sscore output/Simulate_Genotypes/4PopSplit/V4/genos.pop --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules simulatate_phenotype_4PopSplit --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/152.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/152.jobfailed; exit 1)

