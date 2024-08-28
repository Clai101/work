#!/bin/bash

rm hbk_data/*.hist -f
rm log_data/*.log -f
rm index_data/*.index -f
rm -f fpda_pid.*

make

for i in  07 09  11 13 15 17 19 21 23 25 27  31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73

do 

./scripts/exp.sh $i &

done

wait

echo "Done!"
