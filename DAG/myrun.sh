#!/bin/bash
JOBNAME=JOBA
DURATION=30
echo Job : $JOBNAME start
# declare a name for this job to be sample_job
#PBS -N JOBNAME
#PBS -l nodes=1:ppn=2,pvmem=2000mb
#PBS -l walltime=1:00:00
# combine PBS standard output and error files
##PBS -j oe
##PBS -k eo

#How many procs do I have?
#NN=`cat $PBS_NODEFILE | wc -l`
#echo "Processors received = "$NN
#echo "script running on host `hostname`"

#cd into the directory where I typed qsub
#cd $PBS_O_WORKDIR
#echo "PBS_NODEFILE"
#cat $PBS_NODEFILE

#Type in commands to run. Replace a.out with the program name
#to run.

//opt/openmpi/bin/mpiexec -n 2 //share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python helloworld.py  $DURATION