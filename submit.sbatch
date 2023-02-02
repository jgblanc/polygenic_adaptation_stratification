#!/bin/bash
#SBATCH --job-name=4PopSplitB
#SBATCH --output=logs/4PopSplitB.out
#SBATCH --error=logs/4PopSplitB.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=20
#SBATCH --mem-per-cpu=2000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
#snakemake --cores all --rerun-incomplete  
snakemake -j6 --rerun-incomplete 
