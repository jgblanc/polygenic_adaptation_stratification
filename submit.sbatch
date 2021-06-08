#!/bin/bash
#SBATCH --job-name=F3-en
#SBATCH --output=logs/F3-en.out
#SBATCH --error=logs/F3-en.err
#SBATCH --time=30:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --mem=8G

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
snakemake --cores 
