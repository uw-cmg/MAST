#!/bin/bash
#PBS -N ?mast_name?
#PBS -q ?mast_queue?
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
LD_LIBRARY_PATH=//share/apps/mvapich2/tam_mvapich2-1.9a2/usr/local/lib:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH
?mast_exec?
