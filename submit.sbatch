#!/bin/bash
#SBATCH --job-name=CalcQ
#SBATCH --output=logs/CalcQ.out
#SBATCH --error=logs/CalcQ.err
#SBATCH --time=07:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mem-per-cpu=4000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete  
