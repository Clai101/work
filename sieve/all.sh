#!/bin/bash

rm hbk_mc/*.hist -f
rm hbk_data/*.hist -f
rm log_mc/*.log -f
rm log_data/*.log -f
rm index_mc/*.index -f
rm index_data/*.index -f
rm -f fpda_pid.*

> scripts/started.txt

make

for i in  07 09  11 13 15 17 19 21 23 25 27  31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73

do 

./scripts/exp.sh ${i}  \>/dev/null "&" 

done

for((j=0;j<29;j++))

do 

 bsub -ql   ./scripts/mc.sh 7 $j

done

for((j=0;j<13;j++))

do

 bsub -ql   ./scripts/mc.sh 9 $j

done

for((j=0;j<14;j++))

do

 bsub -ql   ./scripts/mc.sh 11 $j

done

for((j=0;j<17;j++))

do

 bsub -ql   ./scripts/mc.sh 13 $j

done

for((j=0;j<15;j++))

do

 bsub -ql  ./scripts/mc.sh 15 $j

done

for((j=0;j<11;j++))

do

 bsub -ql   ./scripts/mc.sh 17 $j

done

for((j=0;j<18;j++))

do

 bsub -ql   ./scripts/mc.sh 19 $j

done

for((j=0;j<4;j++))

do

bsub -ql   ./scripts/mc.sh 21 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./scripts/mc.sh 23 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./scripts/mc.sh 25 $j

done

for((j=0;j<17;j++))

do

bsub -ql   ./scripts/mc.sh 27 $j

done

for((j=0;j<18;j++))

do

bsub -ql   ./scripts/mc.sh 31 $j

done

for((j=0;j<9;j++))

do

bsub -ql   ./scripts/mc.sh 33 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./scripts/mc.sh 35 $j

done

for((j=0;j<20;j++))

do

bsub -ql   ./scripts/mc.sh 37 $j

done

for((j=0;j<14;j++))

do

bsub -ql   ./scripts/mc.sh 39 $j

done


for((j=0;j<10;j++))

do

bsub -ql   ./scripts/mc.sh 43 $j

done

for((j=0;j<10;j++))

do

bsub -ql   ./scripts/mc.sh 53 $j

done

for((j=0;j<7;j++))

do

bsub -ql   ./scripts/mc.sh 67 $j

done

for((j=0;j<14;j++))

do

bsub -ql   ./scripts/mc.sh 69 $j

done

for((j=0;j<22;j++))

do

bsub -ql   ./scripts/mc.sh 71 $j

done


