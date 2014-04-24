#!/bin/bash
#PBS -N ?mast_name?
#PBS -l nodes=?mast_nodes?:ppn=?mast_ppn?,pvmem=?mast_memory?mb
#PBS -l walltime=?mast_walltime?
##export all environment variables. 
#PBS -V
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
?mast_exec?
