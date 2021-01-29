#!/bin/sh
# properties = {"type": "single", "rule": "all", "local": true, "input": ["output/PRS/4PopSplit/V2/C1/genos-test_common.c.p.sscore", "output/PRS/4PopSplit/V2/C2/genos-test_common.c.p.sscore", "output/PRS/4PopSplit/V2/C1/genos-test_common.true.sscore", "output/PRS/4PopSplit/V2/C2/genos-test_common.true.sscore"], "output": [], "wildcards": {}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 0, "cluster": {"mem": 16000, "partition": "broadwl", "ntasks": 1, "tasks": 1, "mem-per-cpu": 2000, "output": "logs/all..out", "error": "logs/all..err", "job-name": "all."}}
 cd /project2/jjberg/jgblanc/polygenic_adaptation_stratification && \
/software/python-3.7.0-el7-x86_64/bin/python \
-m snakemake all --snakefile /project2/jjberg/jgblanc/polygenic_adaptation_stratification/snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs output/PRS/4PopSplit/V2/C1/genos-test_common.c.p.sscore output/PRS/4PopSplit/V2/C2/genos-test_common.c.p.sscore output/PRS/4PopSplit/V2/C1/genos-test_common.true.sscore output/PRS/4PopSplit/V2/C2/genos-test_common.true.sscore --latency-wait 5 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
  -p --allowed-rules all --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/0.jobfinished || (touch /project2/jjberg/jgblanc/polygenic_adaptation_stratification/.snakemake/tmp.zd_ltufs/0.jobfailed; exit 1)

