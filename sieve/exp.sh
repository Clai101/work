#!/bin/bash
if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi


for ((i=0; i<30; i++)) do

for ((j=0; j<10; j++)) do
if [ `./check_statistics.sh  $1 $i $j` -ne "0" ] 
then 
bsub -ql ./data.sh ${1} ${i} ${j}
fi
done
done

