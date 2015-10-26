import os
try:
    from ase import Atom, Atoms
    from ase.calculators.vasp import Vasp
    from ase.optimize import BFGS
    from ase.units import GPa
    from ase.calculators.neighborlist import NeighborList
except ImportError:
    print "NOTE: ASE is not installed. To use Structopt eval_energy.py, ASE must be installed."
from MAST.structopt.inp_out.write_xyz import write_xyz
from MAST.structopt.tools.setup_calculator import setup_calculator
from MAST.structopt.tools.find_defects import find_defects
from MAST.structopt.tools.check_cell_type import check_cell_type
from MAST.structopt.fingerprinting import get_fingerprint
from MAST.structopt.tools.lammps import LAMMPS
import numpy
import math
try:
    from mpi4py import MPI
except ImportError:
    print "NOTE: mpi4py is not installed. To use certain features in Structopt eval_energy.py, mpi4py must be installed."
import logging
import pdb
import shutil
import time
import scipy
import random

def eval_energy(Optimizer, individ):
    """Function to evaluate energy of an individual
    Inputs:
        input = [Optimizer class object with parameters, Individual class structure to be evaluated]
    Outputs:
        energy, bul, individ, signal
        energy = energy of Individual evaluated
        bul = bulk structure of Individual if simulation structure is Defect
        individ = Individual class structure evaluated
        signal = string of information about evaluation
    """
    if 'stem' in Optimizer.fitness_scheme:
        return structopt.tools.eval_energy_stem(Optimizer, individ)
    else:
        return structopt.tools.eval_energy_non_stem(Optimizer, individ)
