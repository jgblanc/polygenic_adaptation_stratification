#!/bin/sh

a=1

while [ $a -lt 1001 ]
do
   echo $a
   Rscript generate_prs_files.R "A${a}" 
   a=`expr $a + 1`
done

