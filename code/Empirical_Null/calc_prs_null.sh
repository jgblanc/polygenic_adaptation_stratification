#!/bin/bash

pfile=$1
freq=$2
betas_pre=$3
out_pre=$4
num=$5

i=1
while [ "$i" -le "$num" ]; do

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.c.betas" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.c"  \
   --score-col-nums 3,4

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.c.p.betas" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.c.p"  \
   --score-col-nums 3,4

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.nc.betas" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.nc"  \
   --score-col-nums 3,4

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.c.betas-Tm" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.c-Tm"  \
   --score-col-nums 3,4

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.c.p.betas-Tm" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.c.p-Tm"  \
   --score-col-nums 3,4

  plink2 \
   --pfile $pfile \
   --read-freq $freq \
   --score "${betas_pre}${i}.nc.betas-Tm" cols=dosagesum,scoresums \
   --out "${out_pre}${i}.nc-Tm"  \
   --score-col-nums 3,4

  i=$(($i + 1))
done