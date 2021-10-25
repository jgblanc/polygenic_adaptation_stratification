#!/bin/bash
#SBATCH --job-name=4PopSplit_T100h23
#SBATCH --output=logs/4PopSplit_T100h23.out
#SBATCH --error=logs/4PopSplit_T100h23.err
#SBATCH --time=05:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=6GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake -j 20 --rerun-incomplete    
#sh copy_dir.sh
