#!/bin/bash
#SBATCH --job-name=4PopSplit_K100
#SBATCH --output=logs/4PopSplit_K100.out
#SBATCH --error=logs/4PopSplit_K100.err
#SBATCH --time=02:00:00
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
snakemake --cores 4
