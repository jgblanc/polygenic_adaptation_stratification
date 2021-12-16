#!/bin/bash
#SBATCH --job-name=SimpleGrid_PopID
#SBATCH --output=logs/SimpleGrid_PopID.out
#SBATCH --error=logs/SimpleGrid_PopID.err
#SBATCH --time=04:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=5
#SBATCH --ntasks-per-node=4
#SBATCH --mem-per-cpu=6GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake -j 20 --rerun-incomplete    
#sh copy_dir.sh
