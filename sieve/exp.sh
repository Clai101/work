#!/bin/bash
if test $1
then
    exp=$1
else
    echo exp number is missing
fi


for ((i=0; i<30; i++)) do

    for ((j=0; j<10; j++)) do
        if [ `./check_statistics.sh  $1 $i $j` -ne "0" ] 
        then 
            echo "data_${1}_${i}_${j}" >> started.txt
            bsub -qs ./data.sh ${1} ${i} ${j}
        fi
    done
done
