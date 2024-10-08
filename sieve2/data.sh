#!/bin/bash

if test $1
then
    exp=$1
else
    echo exp number is missing
    exit 1
fi

if test $2
then
    run=$2
    else
    echo run number is missing
    exit 1
fi

source /sw/belle/local/etc/bashrc_general

export USE_GRAND_REPROCESS_DATA=1
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so

#export BELLE_MESSAGE_LEVEL=DDEBUG
export BELLE_MESSAGE_LEVEL=INFO
#### case B

unset BELLE_USE_TMP

filelist=`find /gpfs/home/belle2/matrk/work/sieve/index_data -name $exp.$run\*.index -size +0c|sort`

outfile=hbk_data/$exp.$run.hist
logfile=log_data/$exp.$run.log

(
cat <<EOF

module register fix_mdst User_reco user_index
path add_module main fix_mdst User_reco user_index
path add_condition main <:0:KILL

initialize
output open       index_data/${exp}.${run}.index
histogram define $outfile

EOF

for file in ${filelist}
do
    echo "process_event ${file} " 
done

echo "terminate"
) |basf> $logfile 2>&1
