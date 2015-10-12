#!/bin/sh
mast_head_dir="//home/mayeshiba/testing_MAST"
script_head_dir="//mnt/gluster/mayeshiba/MAST/MAST/submit/platforms/chtc_dagman"
output_file="PRE_script_output"

date >> $output_file
export MAST_SCRATCH=$mast_head_dir"/SCRATCH"
export MAST_ARCHIVE=$mast_head_dir"/ARCHIVE"
export MAST_CONTROL=$mast_head_dir"/CONTROL"
export MAST_PLATFORM="chtc_dagman"
python $script_head_dir"/mast_do_setup_python.py" $MAST_SCRATCH"/"$1 $2 1>> $output_file 2>> $output_file

echo $2
if [[ $2 == *"neb_"* ]]
then 
    echo "NEB detected; switching submission script" >> $output_file
    cp $script_head_dir"/submit_neb.sh" $MAST_SCRATCH"/"$1"/"$2"/submit.sh"
else echo "not neb" >> $output_file
fi
