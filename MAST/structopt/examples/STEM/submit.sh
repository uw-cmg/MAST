#!/bin/bash
#PBS -N test_Au
#PBS -q morgan2
#PBS -l nodes=1:ppn=12,pvmem=1000mb
#PBS -l walltime=24:00:00
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
LAMMPS_COMMAND="~/lmp_serial"
export LAMMPS_COMMAND
/share/apps/mvapich2/tam_mvapich2-1.9a2/usr/local/bin/mpirun /share/apps/EPD_64bit/epd_free-7.3-2-rh5-x86_64/bin/python STEM_run.py > output.txt









































































