#!/bin/bash
#SBATCH --job-name=snkmk-big
#SBATCH --output=logs/snkmk-big.out
#SBATCH --error=logs/snkmk-big.err
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
snakemake --cores 8 aggregate_genotypes
