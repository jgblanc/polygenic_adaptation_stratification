#!/bin/bash
#SBATCH --job-name=SimpleGrid_T100_stratD
#SBATCH --output=logs/SimpleGrid_T100_stratD.out
#SBATCH --error=logs/SimpleGrid_T100_stratD.err
#SBATCH --time=20:00:00
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
