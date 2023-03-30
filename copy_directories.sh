#!/bin/bash
for i in {3..100}
do
   echo "Print $i times"
   cp -R output/Simulate_Genotypes/4PopSplit/C${i}/*  output/Simulate_Genotypes/4PopSplit/B${i}/
done
