#!/bin/bash

exp=$1
run=$2
run1=$3                                                                                                                                                                                               
source /sw/belle/local/etc/bashrc_general
(
check_process_url "http://bweb3/montecarlo.php?bl=caseB&ty=evtgen-charm&ex=${exp}&rs=${run}${run1}0&re=${run}${run1}9&dv=zfserv"
) > mcscripts/${exp}.${run}.${run1}.script 2>&1
