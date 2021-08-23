#!/bin/bash
#SBATCH --job-name=4PopSplit_F100_rescale
#SBATCH --output=logs/4PopSplit_F100_rescale.out
#SBATCH --error=logs/4PopSplit_F100_rescale.err
#SBATCH --time=24:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=10
#SBATCH --mem-per-cpu=4GB

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 10 --rerun-incomplete    
