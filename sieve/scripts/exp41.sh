#!/bin/bash

for ((i=0; i<13; i++)) do

for ((j=0; j<10; j++)) do
if [ `../check_statistics.sh  41 $i $j` -ne "0" ] 
then

bsub -ql ../data.sh 41 ${i} ${j}
fi
done
done
 
