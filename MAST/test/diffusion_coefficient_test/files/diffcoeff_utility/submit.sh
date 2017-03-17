#!/bin/bash
#PBS -N diffcoeff_utility
#PBS -l nodes=1:ppn=1,pvmem=1000mb
#PBS -l walltime=4:00:00
##export all environment variables. 
#PBS -V
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
python /home/tam/tammast/MAST/utility/DiffusionCoefficient.py -i diffcoeff_input.txt
