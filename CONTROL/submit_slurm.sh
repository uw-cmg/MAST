#!/bin/bash

#SBATCH -t 1:00:00
#SBATCH -J mastmon
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --ntasks-per-node=1
python runmast.py
