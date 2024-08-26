#!/bin/bash

rm hbk_mc/*.hist -f
rm hbk_data/*.hist -f
rm log_mc/*.log -f
rm log_data/*.log -f
rm index_mc/*.index -f
rm index_data/*.index -f
rm -f fpda_pid.*

> started.txt

make

for i in  07 09  11 13 15 17 19 21 23 25 27 31 33 35 37 39 41 43 45 47 49 51 53 55 61 63 65 67 69 71 73

do 

./exp.sh ${i}

done


for((j=0;j<29;j++))

do 

    echo "mc_7_${j}" >> started.txt
    bsub -ql   ./mc.sh 7 $j

done

for((j=0;j<13;j++))

do

    echo "mc_9_${j}" >> started.txt
    bsub -ql   ./mc.sh 9 $j

done

for((j=0;j<14;j++))

do

    echo "mc_11_${j}" >> started.txt
    bsub -ql   ./mc.sh 11 $j

done

for((j=0;j<17;j++))

do

    echo "mc_13_${j}" >> started.txt
    bsub -ql   ./mc.sh 13 $j

done

for((j=0;j<15;j++))

do

    echo "mc_15_${j}" >> started.txt
 bsub -ql  ./mc.sh 15 $j

done

for((j=0;j<11;j++))

do

    echo "mc_17_${j}" >> started.txt
    bsub -ql   ./mc.sh 17 $j

done

for((j=0;j<18;j++))

do

    echo "mc_19_${j}" >> started.txt
    bsub -ql   ./mc.sh 19 $j

done

for((j=0;j<4;j++))

do

    echo "mc_21_${j}" >> started.txt
    bsub -ql   ./mc.sh 21 $j

done

for((j=0;j<7;j++))

do

    echo "mc_23_${j}" >> started.txt
    bsub -ql   ./mc.sh 23 $j

done

for((j=0;j<22;j++))

do

    echo "mc_25_${j}" >> started.txt
    bsub -ql   ./mc.sh 25 $j

done

for((j=0;j<17;j++))

do

    echo "mc_27_${j}" >> started.txt
    bsub -ql   ./mc.sh 27 $j

done

for((j=0;j<18;j++))

do

    echo "mc_31_${j}" >> started.txt
    bsub -ql   ./mc.sh 31 $j

done

for((j=0;j<9;j++))

do

    echo "mc_33_${j}" >> started.txt
    bsub -ql   ./mc.sh 33 $j

done

for((j=0;j<7;j++))

do

    echo "mc_35_${j}" >> started.txt
    bsub -ql   ./mc.sh 35 $j

done

for((j=0;j<20;j++))

do

    echo "mc_37_${j}" >> started.txt
    bsub -ql   ./mc.sh 37 $j

done

for((j=0;j<14;j++))

do

    echo "mc_39_${j}" >> started.txt
    bsub -ql   ./mc.sh 39 $j

done


for((j=0;j<10;j++))

do

    echo "mc_43_${j}" >> started.txt
    bsub -ql   ./mc.sh 43 $j

done

for((j=0;j<10;j++))

do

    echo "mc_53_${j}" >> started.txt
    bsub -ql   ./mc.sh 53 $j

done

for((j=0;j<7;j++))

do

    echo "mc_67_${j}" >> started.txt
    bsub -ql   ./mc.sh 67 $j

done

for((j=0;j<14;j++))

do

    echo "mc_69_${j}" >> started.txt
    bsub -ql   ./mc.sh 69 $j

done

for((j=0;j<22;j++))

do

    echo "mc_71_${j}" >> started.txt
    bsub -ql   ./mc.sh 71 $j

done


