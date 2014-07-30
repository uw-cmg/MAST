#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -J mastmon
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -p serial
#SBATCH -o mastmon_submission_output
#SBATCH -e mastmon_submission_errors
module load python
module load vasp
python runmast.py >> $MAST_CONTROL/mastoutput 2> $MAST_CONTROL/errormast
