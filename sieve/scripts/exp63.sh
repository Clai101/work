#!/bin/bash

for ((i=0; i<8; i++)) do

for ((j=0; j<10; j++)) do
if [ `./check_statistics.sh  63 $i $j` -ne "0" ] 
then

bsub -ql ./data.sh 63 ${i} ${j}
fi 
done
done

