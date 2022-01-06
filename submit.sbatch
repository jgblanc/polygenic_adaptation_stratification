#!/bin/bash
#SBATCH --job-name=4PopSplit_A_Genotypes
#SBATCH --output=logs/4PopSplit_A_Genotypes.out
#SBATCH --error=logs/4PopSplit_A_Genotypes.err
#SBATCH --time=24:00:00
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
snakemake -j 20 --rerun-incomplete aggregate_genotypes   
#sh copy_dir.sh
