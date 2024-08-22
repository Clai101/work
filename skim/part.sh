#!/bin/bash

rm hbk_mc/* -f
rm hbk_data/* -f
rm log_mc/* -f
rm log_data/* -f
rm index_mc/* -f
rm index_data/* -f
rm -f fpda_pid.*

make

for i in 7 9 11 13 15 17 19 21 23 25 27 31
do 

bsub -ql   ./mc.sh $i
bsub -ql   ./data.sh $i

done