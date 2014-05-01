#!/bin/bash

#SBATCH -t ?mast_walltime?
#SBATCH -J ?mast_name?
#SBATCH -N ?mast_nodes?
#SBATCH -n ?mast_processors?
#SBATCH --ntasks-per-node=?mast_ppn?
?mast_exec?
