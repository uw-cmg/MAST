#!/bin/bash
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2016-02-08
##############################################################
##############################################################
# Arguments:
# [$1] Full path to example files
# [$2] Name of input file
# [$3] Activation command (optional)
# [$4] Activation environment (optional)
##############################################################
examples_located=$1
which_example=$2

if [ $# -gt 2 ]
then
    echo "Sourcing."
    activate_command=$3
    testing_environment=$4
    source $activate_command $testing_environment
    echo "Source done."
fi

timestamp=`date +%Y%m%d"T"%H%M%S`
total_mast_calls=500
sleep_interval=60

export MAST_PLATFORM="no_queue_system"
mast_test=$HOME"/MAST/workflow_test_"$timestamp
mkdir $mast_test
mkdir $mast_test/SCRATCH
mkdir $mast_test/CONTROL
mkdir $mast_test/ARCHIVE
export MAST_SCRATCH=$mast_test/SCRATCH
export MAST_ARCHIVE=$mast_test/ARCHIVE
export MAST_CONTROL=$mast_test/CONTROL
echo "" >> $MAST_CONTROL/submitlist
echo "" >> $MAST_CONTROL/just_submitted
which mast

mast -i $examples_located"/"$which_example
mydir=`find $MAST_SCRATCH -mindepth 1 -maxdepth 1 -type ;d`

cd $mydir

mct="0"

while [ $mct -lt $total_mast_calls ]
do
    date
    echo "MAST call $mct"
    mast
    if [ -f "SUMMARY.txt" ]
    then
        echo "SUMMARY file found. Workflow completed."
        break
    fi
    sleep $sleep_interval
    mct=$[$mct+1]
done



