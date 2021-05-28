#!/bin/bash

prs_pre=$1
lambda_T=$2
Va=$3
Va_Tm=$4
true_file=$5
Tvec=$6
out_pre=$7
num=$8

i=1
while [ "$i" -le "$num" ]; do

   Rscript code/Empirical_Null/calc_Qx.R \
   "${prs_pre}${i}.c.sscore" \
   "${prs_pre}${i}.c.p.sscore" \
   "${prs_pre}${i}.nc.sscore" \
   "${prs_pre}${i}.c-Tm.sscore" \
   "${prs_pre}${i}.c.p-Tm.sscore" \
   "${prs_pre}${i}.nc-Tm.sscore" \
   "${lambda_T}" \
   "${Va}" \
   "${Va_Tm}" \
   "${true_file}" \
   "${Tvec}" \
   "${out_pre}${i}.txt"

  i=$(($i + 1))
done