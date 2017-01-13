#!/bin/bash
#PBS -N neb_ring1-ring2_q=p0_opt1
#PBS -q morgan3
#PBS -l nodes=1:ppn=32,pvmem=1000mb
#PBS -l walltime=96:00:00
##export all environment variables. 
#PBS -V
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
//opt/mpiexec/bin/mpiexec //share/apps/bin/vasp5.2_CNEB
