#!/bin/bash

rm hbk_data/*.hist -f &
rm log_data/*.log -f &
rm index_data/*.index -f &

wait


for i in  07 09  11 13 15
do 

./exp.sh $i data &

done

wait

echo "Done!"
