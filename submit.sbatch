#!/bin/bash
#SBATCH --job-name=snkmk-en1
#SBATCH --output=logs/snkmk-en1.out
#SBATCH --error=logs/snkmk-en1.err
#SBATCH --time=08:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --mem=12G

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
snakemake --cores 8 
