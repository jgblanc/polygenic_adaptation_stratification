#!/bin/bash
#SBATCH --job-name=4PopSplit_Concat
#SBATCH --output=logs/4PopSplit_Concat.out
#SBATCH --error=logs/4PopSplit_Concat.err
#SBATCH --partition=bigmem2
#SBATCH --time=09:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem-per-cpu=16000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
#snakemake --cores all --rerun-incomplete  
snakemake -j3 --rerun-incomplete 
