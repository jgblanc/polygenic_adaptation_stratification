#!/bin/bash
for i in {3..100}
do
   echo "Welcome $i times"
   cp output/Simulate_Genotypes/4PopSplit/C${i}/genos*  output/Simulate_Genotypes/4PopSplit/A${i}/
done
