#!/bin/bash
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2016-02-08
##############################################################
##############################################################
# Arguments:
# [$1] MAST testing tree to use
# [$2] Full path to example files
# [$3] Name of input file
# [$4] Activation command (optional)
# [$5] Activation environment (optional)
##############################################################
mast_test_dir=$1
examples_located=$2
which_example=$3

if [ $# -gt 3 ]
then
    echo "Sourcing."
    activate_command=$4
    testing_environment=$5
    source $activate_command $testing_environment
    echo "Source done."
fi

timestamp=`date +%Y%m%d"T"%H%M%S`
total_mast_calls=500
sleep_interval=60

export MAST_PLATFORM="no_queue_system"
mast_test=$HOME"/MAST/workflow_test_"$timestamp
export MAST_SCRATCH=$mast_test_dir/SCRATCH
export MAST_ARCHIVE=$mast_test_dir/ARCHIVE
export MAST_CONTROL=$mast_test_dir/CONTROL

which mast

mast -i $examples_located"/"$which_example
mydir=`find $MAST_SCRATCH -mindepth 1 -maxdepth 1 -type d`

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
        echo "Workflow at:"
        echo $mydir
        break
    fi
    sleep $sleep_interval
    mct=$[$mct+1]
done



