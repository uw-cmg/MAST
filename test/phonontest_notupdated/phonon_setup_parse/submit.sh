#!/bin/bash
#PBS -N phonontest_phonon_vac1_parse
#PBS -l nodes=1:ppn=1,pvmem=1000mb
#PBS -l walltime=4:00:00
NN=`cat $PBS_NODEFILE | wc -l`
echo "Processors received = "$NN
echo "script running on host `hostname`"
cd $PBS_O_WORKDIR
echo "PBS_NODEFILE"
cat $PBS_NODEFILE
/home/tam/tammast/bin/phon_nosymm
