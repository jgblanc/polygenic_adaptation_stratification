#!/bin/sh

for i in {3..100}
do
   echo "$i"
   cp -r output/Simulate_Genotypes/SimpleGrid/T${i}/C1/ output/Simulate_Genotypes/SimpleGrid/T${i}/C3/
done
