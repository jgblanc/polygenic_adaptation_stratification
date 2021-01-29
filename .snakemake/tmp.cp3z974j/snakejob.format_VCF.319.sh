#!/bin/sh
# properties = {"type": "single", "rule": "format_VCF", "local": false, "input": ["output/Simulate_Genotypes/4PopSplit/V8/genos_7.vcf"], "output": ["output/Simulate_Genotypes/4PopSplit/V8/genos_7.ids.vcf.gz"], "wildcards": {"model": "4PopSplit", "rep": "V8", "chr": "7"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 319, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/format_VCF.chr=7,model=4PopSplit,rep=V8.out", "error": "logs/format_VCF.chr=7,model=4PopSplit,rep=V8.err", "job-name": "format_VCF.chr=7,model=4PopSplit,rep=V8"}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake output/Simulate_Genotypes/4PopSplit/V8/genos_7.ids.vcf.gz --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j output/Simulate_Genotypes/4PopSplit/V8/genos_7.vcf --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules format_VCF --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/319.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.cp3z974j/319.jobfailed; exit 1)

