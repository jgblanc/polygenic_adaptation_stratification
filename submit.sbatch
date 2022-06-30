#!/bin/bash
#SBATCH --job-name=4PopSplit_A_Shift
#SBATCH --output=logs/4PopSplit_A_Shift.out
#SBATCH --error=logs/4PopSplit_A_Shift.err
#SBATCH --time=01:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=14
#SBATCH --mem-per-cpu=4000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete  
