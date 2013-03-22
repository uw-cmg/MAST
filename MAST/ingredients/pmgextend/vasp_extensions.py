from pymatgen.io.vaspio import *
import numpy as np

def get_max_enmax_from_potcar(mypotcar):
    """Get maximum enmax value (float) from Potcar (combined list)"""
    enmax_list=list()
    potcarct=0
    onepotcar=None
    while potcarct < len(mypotcar):
        onepotcar = mypotcar[potcarct] #A PotcarSingle object
        enmax_list.append(onepotcar.enmax)
        potcarct = potcarct + 1
    return max(enmax_list)

def make_one_unfrozen_atom_poscar(myposcar, natom):
    """Use selective dynamics to make a poscar with one unfrozen atom.
        myposcar = Poscar
        natom = the number of the atom to unfreeze
        Returns: Poscar (use write_file function on it).
    """
    mysd=np.zeros([sum(myposcar.natoms),3],bool)
    mysd[natom-1][0]=True #indexing starts at 0
    mysd[natom-1][1]=True
    mysd[natom-1][2]=True
    myposcar.selective_dynamics = mysd
    return myposcar

def make_one_unfrozen_direction_poscar(myposcar, natom, ndir):
    """Use selective dynamics to make a poscar with one unfrozen atom.
        myposcar = Poscar
        natom = the number of the atom to unfreeze
        ndir = the direction to freeze (0, 1, 2 for x, y, z)
        Returns: Poscar (use write_file function on it).
    """
    mysd=np.zeros([sum(myposcar.natoms),3],bool)
    mysd[natom-1][ndir]=True #indexing starts at 0
    myposcar.selective_dynamics = mysd
    return myposcar
