#!/bin/bash

for ((i=0; i<24; i++)) do

for ((j=0; j<10; j++)) do
if [ `../check_statistics.sh  71 $i $j` -ne "0" ] 
then

bsub -q s ../data.sh 71 ${i} ${j} 
fi
done
done

