#!/bin/sh

date >> PRE_script_output
export MAST_PLATFORM="chtc_dagman"
python //home/mayeshiba/DAG_for_MAST/mast_do_setup_python.py $MAST_SCRATCH"/"$1 $2 1>> PRE_script_output 2>> PRE_script_output
