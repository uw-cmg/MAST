#!/bin/bash
##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2016-02-08
##############################################################
##############################################################
# Requirements:
# 1. Home directory access from where the test will be run
# 2. MAST installation
##############################################################
##############################################################
# Arguments:
# [$1] Input file name
##############################################################
######################################
# Modify this section as necessary
######################################
activate_command="$HOME/anaconda/bin/activate"
testing_environment="mast_install_20160115"
examples_located="$HOME/mast_2016_tam/MAST/examples"
testdir="$HOME/mast_2016_tam/MAST/test/workflow_test"
#######

which_example = $1

#create unique MAST tree for testing
timestamp=`date +%Y%m%d"T"%H%M%S`
mast_test_dir=$HOME"/MAST/workflow_test_"$timestamp
mkdir $mast_test_dir
cp -r $testdir/mini_mast_tree/* $mast_test_dir/.
#submit workflow test to queue
shortname=`echo $which_example | awk -F. '{print $1}'`
output="output_"$shortname
submitscript="submit_"$shortname".sh"
bashcommand="bash $testdir/generic_mast_workflow.sh $mast_test_dir $examples_located $which_example $activate_command $testing_environment >> $output"
cp $testdir/submit_stub.sh $testdir/$submitscript
echo $bashcommand >> $testdir/$submitscript
qsub $testdir/$submitscript

echo $mast_test_dir
