#!/bin/sh

for i in {55,90,71,66,96,51,97,78,83,68,74,99,100}
do
   echo "$i"
   bgzip output/Calculate_Tm/SimpleGrid/T${i}/C2/pca.eigenvec.allele
   bgzip output/Calculate_Tm/SimpleGrid/T${i}/C2/projection.sscore
done
