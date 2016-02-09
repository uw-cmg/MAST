#!/bin/bash
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2016-02-08
##############################################################
##############################################################
# This test only works where there is a shared home directory.
# 1) For initial testing (unsure of run),
# use interactive submission to get on a node, e.g.
#     qsub -I -q morganshort
#     qrsh -q morganeth.q
#     srun -n1 -N1 -p int --pty bash
# 2) Otherwise, if fairly sure of test, submit this script to a queue
##############################################################

######################################
# Modify this section as necessary
######################################
activate_command="$HOME/anaconda/bin/activate"
testing_environment="mast_install_20160115"
examples_located="$HOME/mast_2016_tam/MAST/examples"
which_example="simple_optimization.inp"
folder_name_match="OptimizeWorkflowTest*"
#seconds between MAST calls
sleep_interval="30"
total_mast_calls="30"
#######

timestamp=`date +%Y%m%d"T"%H%M%S`
source $activate_command $testing_environment
export MAST_PLATFORM="no_queue_system"
mkdir $HOME"/MAST/scratch_test_"$timestamp
export MAST_SCRATCH=$HOME"/MAST/scratch_test_"$timestamp
which mast

mast -i $examples_located"/"$which_example

cd $MAST_SCRATCH"/"$folder_name_match

mct="0"

while [ $mct -lt $total_mast_calls ]
do
    date
    echo "MAST call $mct"
    mast
    if [ -f "SUMMARY.txt" ]
    then
        break
    fi
    sleep $sleep_interval
    mct=$[$mct+1]
done



