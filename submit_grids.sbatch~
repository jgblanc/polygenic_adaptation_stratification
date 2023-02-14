#!/bin/bash
#SBATCH --job-name=4PopSplit_Reps3
#SBATCH --output=logs/4PopSplit_Rep3.out
#SBATCH --error=logs/4PopSplit_Rep3.err
#SBATCH --time=09:00:00
#SBATCH --partition=broadwl
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=2000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
snakemake --cores all --rerun-incomplete -s snakefile_reps3  
#snakemake -j1 --rerun-incomplete  -s snakefile_reps2 
