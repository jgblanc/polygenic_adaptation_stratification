#!/bin/bash
#SBATCH --job-name=F102
#SBATCH --output=logs/F102.out
#SBATCH --error=logs/F102.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --mem=4G

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
snakemake --cores 8 --rerun-incomplete
