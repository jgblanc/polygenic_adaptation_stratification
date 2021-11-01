#!/bin/bash
#SBATCH --job-name=SimpleGrid_PS_tv
#SBATCH --output=logs/SimpleGrid_PS_tv.out
#SBATCH --error=logs/SimpleGrid_PS_tv.err
#SBATCH --time=08:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake -j 20 --rerun-incomplete    
#sh copy_dir.sh
