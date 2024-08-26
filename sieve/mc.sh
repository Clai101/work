#!/bin/sh

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

#####################################
source /sw/belle/local/etc/bashrc_general


export USE_GRAND_REPROCESS_DATA=1
export BASF_USER_IF=basfsh.so
export BASF_USER_INIT=user_init.so



#export BELLE_MESSAGE_LEVEL=DDEBUG
export BELLE_MESSAGE_LEVEL=INFO
#### case B

unset BELLE_USE_TMP

(
cat <<EOF

module register fix_mdst User_reco user_index
path add_module main fix_mdst  User_reco user_index
path add_condition main <:0:KILL

initialize
output open       index_mc/${exp}.${run}.index
histogram define  hbk_mc/${exp}.${run}.hist

process_url "http://bweb3/montecarlo.php?bl=caseB&ty=evtgen-charm&ex=${exp}&rs=${run}00&re=${run}99&dv=zfserv"

EOF

echo terminate

) | basf >  log_mc/${exp}.${run}.log 2>&1

