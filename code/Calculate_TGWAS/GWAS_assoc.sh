#!/bin/bash
pfile_path=$1
pheno_path=$2
outfile=$3

plink2 \
  --pfile $pfile_path \
  --glm allow-no-covars omit-ref \
  --pheno $pheno_path \
  --pheno-name pheno \
  --geno-counts \
  --out $outfile
