#!/bin/bash

rm hbk_data/*.hist -f &
rm log_data/*.log -f &
rm index_data/*.index -f &

wait


for i in  07 09  11 13 15 17 19 21 23 25 27  31 33 35
do 

./exp.sh $i data &

done

wait

echo "Done!"
