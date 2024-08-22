#!/bin/sh

for i in 7 9 11 13 15 17 19 21 23 25 27 31 33 37 39 43 53 67 69 71
do 

bsub -qs   ./index.sh $i

done

