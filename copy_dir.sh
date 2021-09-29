#!/bin/sh

for i in {51..100}
do
   echo "$i"
   bgzip output/Calculate_Tm/SimpleGrid/T${i}/C2/pca.eigenvec.allele
   bgzip output/Calculate_Tm/SimpleGrid/T${i}/C2/projection.sscore
done
