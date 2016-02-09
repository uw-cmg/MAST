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

######################################
# Modify this section as necessary
######################################
activate_command="$HOME/anaconda/bin/activate"
testing_environment="mast_install_20160115"
examples_located="$HOME/mast_2016_tam/MAST/examples"
testdir="$HOME/mast_2016_tam/MAST/test/workflow_tests"
#######


for which_example in "simple_optimization.inp"
do
    shortname=`echo $which_example | awk -F. '{print $1}'`
    output="output_"shortname
    submitscript="submit_"shortname".sh"
    bashcommand="bash $testdir/generic_mast_workflow.sh $examples_located $which_example $activate_command $testing_environment >> $output"
    cp $testdir/submit_stub.sh $testdir/$submitscript
    echo $bashcommand >> $testdir/$submitscript
    qsub $testdir/$submitscript
done
