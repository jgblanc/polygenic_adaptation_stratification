#!/bin/sh

a=1

while [ $a -lt 1001 ]
do
   mkdir -p "../../output/PRS/4PopSplit/A${a}"
   mkdir -p "../../output/PRS/4PopSplit/A${a}/C1"
   mkdir -p "../../output/PRS/4PopSplit/A${a}/C1/h2-0"
   mkdir -p "../../output/PRS/4PopSplit/A${a}/C1/h2-0/env-0"
   a=`expr $a + 1`
done


a=1

while [ $a -lt 1001 ]
do
   mkdir -p "../../output/Calculate_Tm/4PopSplit/A${a}"
   mkdir -p "../../output/Calculate_Tm/4PopSplit/A${a}/C1"
   a=`expr $a + 1`
done

a=1

while [ $a -lt 1001 ]
do
   mkdir -p "../../output/Simulate_Genotypes/4PopSplit/A${a}"
   mkdir -p "../../output/Simulate_Genotypes/4PopSplit/A${a}/C1"
   a=`expr $a + 1`
done

