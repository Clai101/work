#!/bin/bash

for ((i=0; i<12; i++)) do

for ((j=0; j<10; j++)) do
if [ `../check_statistics.sh  43 $i $j` -ne "0" ]
then

bsub -q s ../data.sh 43 ${i} ${j} 
fi
done
done
