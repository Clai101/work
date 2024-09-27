#!/bin/bash

rm hbk_mc/*.hist -f &
rm log_mc/*.log -f &
rm index_mc/*.index -f &

wait


for i in  7 9  11 13 15 17 19 21 23 25 27  31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73

do 

./exp.sh $i mc &

done

wait

echo "Done!"
