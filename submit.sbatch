#!/bin/bash
#SBATCH --job-name=SimpleGrid_T100_envs
#SBATCH --output=logs/SimpleGrid_T100_envs.out
#SBATCH --error=logs/SimpleGrid_T100_envs.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --mem-per-cpu=4GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 4 --rerun-incomplete  
