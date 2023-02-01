#!/bin/bash
pfile_path=$1
beta_path=$2
outfile=$3

plink2 \
  --pfile $pfile_path \
  --score $beta_path center header-read cols=dosagesum,scoresums list-variants \
  --out $outfile
