#!/bin/bash
pfile_path=$1
pheno_path=$2
outfile=$3

plink2 \
  --pfile $pfile_path \
  --glm omit-ref allow-no-covars \
  --pheno $pheno_path \
  --pheno-name Tvec \
  --geno-counts \
  --out $outfile
