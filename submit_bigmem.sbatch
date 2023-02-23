#!/bin/bash
#SBATCH --job-name=2PopSplit_Plink
#SBATCH --output=logs/2PopSplit_Plink.out
#SBATCH --error=logs/2PopSplit_Plink.err
#SBATCH --partition=bigmem2
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=16000

module load python/cpython-3.7.0
module load R
module load htslib/1.4.1
module load samtools
module load bcftools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile
#snakemake --cores all --rerun-incomplete  
snakemake -j8 --rerun-incomplete -s snakefile_2Pop
