#!/bin/bash

for ((i=0; i<14; i++)) do

for ((j=0; j<10; j++)) do
if [ `../check_statistics.sh  39 $i $j` -ne "0" ] 
then

bsub -q s ../data.sh 39 ${i} ${j} 
fi 
done
done
