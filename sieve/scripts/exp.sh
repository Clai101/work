#!/bin/bash
if test $exp
then
    exp=$exp
else
    echo exp number is missing
    exit 1
fi

for ((i=0; i<30; i++)) do

for ((j=0; j<10; j++)) do
if [ `./check_statistics.sh  ${exp} $i $j` -ne "0" ] 
then 
bsub -ql ./data.sh ${exp} ${i} ${j}
fi
done
done

