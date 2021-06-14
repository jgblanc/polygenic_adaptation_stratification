#!/bin/bash
#SBATCH --job-name=SimpleGrid_F100
#SBATCH --output=logs/SimpleGrid_F100.out
#SBATCH --error=logs/SimpleGrid_F100.err
#SBATCH --time=05:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=8
#SBATCH --mem=4G

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 8
