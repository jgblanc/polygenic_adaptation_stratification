#!/bin/sh

for i in {78..100}
do
   echo "$i"
   bgzip T${i}/C3/pca.eigenvec.allele
   bgzip T${i}/C3/projection.sscore
   
done
