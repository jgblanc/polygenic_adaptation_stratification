#!/bin/sh

a=1

while [ $a -lt 51 ]
do
   echo $a
   Rscript generate_prs_files_ms.R "M${a}"  
   a=`expr $a + 1`
done

