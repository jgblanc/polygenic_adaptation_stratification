#!/bin/sh

for i in {3..100}
do
   echo "$i"
   cp -r output/Simulate_Genotypes/4PopSplit/S${i}/ output/Simulate_Genotypes/4PopSplit/F${i}/
   rm output/Simulate_Genotypes/4PopSplit/F${i}/C*/*_common*
done
