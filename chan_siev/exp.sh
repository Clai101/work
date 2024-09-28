#!/bin/bash
0;95;0c
if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi
if test $2
then
    type=$2
else
    echo type number is missing
    exit 1
fi



for ((i=0; i<30; i++)) do

    for ((j=0; j<10; j++)) do
            bsub -qs /sw/belle/local/bin/centos7-exec ./${type}.sh ${1} ${i} ${j}
    done
done
