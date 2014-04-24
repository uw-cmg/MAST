#!/bin/bash
#PBS -N mastmon
#PBS -l nodes=1:ppn=1,pvmem=2000mb
#PBS -q default
#PBS -l walltime=4:00:00
##export all environment variables. 
#PBS -V
#PBS -o mastmon_submission_output
#PBS -e mastmon_submission_errors
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
python runmast.py >> $MAST_CONTROL/mastoutput 2> $MAST_CONTROL/errormast
