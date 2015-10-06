#!/bin/sh

date >> POST_script_output
export MAST_SCRATCH="//home/mayeshiba/testing_MAST/SCRATCH"
export MAST_ARCHIVE="//home/mayeshiba/testing_MAST/ARCHIVE"
export MAST_CONTROL="//home/mayeshiba/testing_MAST/CONTROL"
export MAST_PLATFORM="chtc_dagman"
python //home/mayeshiba/DAG_for_MAST/mast_check_is_complete_python.py $MAST_SCRATCH"/"$1 $2 1>> POST_script_output 2>> POST_script_output
