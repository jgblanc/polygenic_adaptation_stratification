#!/bin/bash
#SBATCH --job-name=SimpleGrid_GWASPC
#SBATCH --output=logs/SimpleGrid_GWASPC.out
#SBATCH --error=logs/SimpleGrid_GWASPC.err
#SBATCH --time=12:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=10GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake -j 2 --rerun-incomplete    
#sh copy_dir.sh
