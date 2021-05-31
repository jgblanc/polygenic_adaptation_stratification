#!/bin/bash
#SBATCH --job-name=snkmk-en4
#SBATCH --output=logs/snkmk-en4.out
#SBATCH --error=logs/snkmk-en4.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --mem-per-cpu=2000
#SBATCH --ntasks-per-node=14

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
snakemake --cores 56 --rerun-incomplete
