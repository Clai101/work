#!/bin/sh
rm hbk_mc/* -f
rm hbk_data/* -f
rm log_mc/* -f
rm log_data/* -f
rm index_mc/* -f
rm index_data/* -f
rm -f fpda_pid.*

make

for((j=0;j<29;j++))

do 

 bsub -ql   ./mc.sh 7 $j
 bsub -ql   ./data.sh 7 $j

done

for((j=0;j<13;j++))

do

 bsub -ql   ./mc.sh 9 $j
 bsub -ql   ./data.sh 9 $j

done

for((j=0;j<14;j++))

do

 bsub -ql   ./mc.sh 11 $j
 bsub -ql   ./data.sh 11 $j

done

for((j=0;j<17;j++))

do

 bsub -ql   ./mc.sh 13 $j
 bsub -ql   ./data.sh 13 $j

done

for((j=0;j<15;j++))

do

 bsub -ql  ./mc.sh 15 $j
 bsub -ql  ./data.sh 15 $j

done

for((j=0;j<11;j++))

do

 bsub -ql   ./mc.sh 17 $j
 bsub -ql   ./data.sh 17 $j

done

for((j=0;j<18;j++))

do

 bsub -ql   ./mc.sh 19 $j
 bsub -ql   ./data.sh 19 $j

done

for((j=0;j<4;j++))

do

bsub -ql   ./mc.sh 21 $j
bsub -ql   ./data.sh 21 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./mc.sh 23 $j
bsub -ql   ./data.sh 23 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./mc.sh 25 $j
bsub -ql   ./data.sh 25 $j

done

for((j=0;j<17;j++))

do

bsub -ql   ./mc.sh 27 $j
bsub -ql   ./data.sh 27 $j

done

for((j=0;j<18;j++))

do

bsub -ql   ./mc.sh 31 $j
bsub -ql   ./data.sh 31 $j

done

