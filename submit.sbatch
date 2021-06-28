#!/bin/bash
#SBATCH --job-name=4PopSplit_T100_envs2
#SBATCH --output=logs/4PopSplit_T100_envs2.out
#SBATCH --error=logs/4PopSplit_T100_envs2.err
#SBATCH --time=02:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=5
#SBATCH --mem=5G

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores 5 --rerun-incomplete 
