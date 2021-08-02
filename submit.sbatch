#!/bin/bash
#SBATCH --job-name=SimpleGrid_T100_C2_2
#SBATCH --output=logs/SimpleGrid_T100_C2_2.out
#SBATCH --error=logs/SimpleGrid_T100_C2_2.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=6
#SBATCH --mem-per-cpu=12GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 6 --rerun-incomplete  
