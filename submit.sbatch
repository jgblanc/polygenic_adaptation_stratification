#!/bin/bash
#SBATCH --job-name=SimpleGrid_T100_3Pop
#SBATCH --output=logs/SimpleGrid_T100_3Pop.out
#SBATCH --error=logs/SimpleGrid_T100_3Pop.err
#SBATCH --time=36:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=10
#SBATCH --mem-per-cpu=6GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 10 --rerun-incomplete    
