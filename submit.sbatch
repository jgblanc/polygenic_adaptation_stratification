#!/bin/bash
#SBATCH --job-name=4PopSplit_diff
#SBATCH --output=logs/4PopSplit_diff.out
#SBATCH --error=logs/4PopSplit_diff.err
#SBATCH --time=12:00:00
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
