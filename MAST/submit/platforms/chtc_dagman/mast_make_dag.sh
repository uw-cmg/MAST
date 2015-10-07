#!/bin/sh
mast_head_dir="//home/mayeshiba/testing_MAST"
script_head_dir="//mnt/gluster/mayeshiba/MAST/MAST/submit/platforms/chtc_dagman"
output_file="make_dag_output"

date >> $output_file
export MAST_SCRATCH=$mast_head_dir"/SCRATCH"
export MAST_ARCHIVE=$mast_head_dir"/ARCHIVE"
export MAST_CONTROL=$mast_head_dir"/CONTROL"
export MAST_PLATFORM="chtc_dagman"
python $script_head_dir"/make_my_dag.py" $1 $script_head_dir 1>> $output_file 2>> $output_file
