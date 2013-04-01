from pymatgen.io.vaspio import *
import numpy as np
from MAST.utility.dirutil import *
from MAST.utility import MASTError

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

def combine_dynmats(myposcar, mydir):
    """Combine DYNMATs into one hessian.
        myposcar = Poscar
        mydir = top directory for DYNMAT files
    """
    natoms = sum(myposcar.natoms)
    myhess=np.zeros([natoms*3, natoms*3])
    #arrange as x1, y1, z1, x2, y2, z2, etc.
    dynmatlist = walkfiles(mydir, 1, 5, "*DYNMAT*")
    if len(dynmatlist) == 0:
        raise MASTError("pmgextend combine_dynmats", "No DYNMATs found under " + mydir)
    opendyn=""
    dynlines=[]
    datoms=0
    dmats=0
    mct=0
    for onedynmat in dynmatlist:
        opendyn = open(onedynmat,'rb')
        dynlines=opendyn.readlines()
        opendyn.close()
        datoms = int(dynlines[0].split()[1])
        dmats = int(dynlines[0].split()[2])
        mycount=2 #starting line
        while mycount < len(dynlines):
            littlemat=[]
            topatom = int(dynlines[mycount].split()[0])
            whichdir = int(dynlines[mycount].split()[1])
            littlemat = dynlines[mycount+1:mycount+datoms]
            act = 0
            while act < datoms:
                dactx=int(littlemat[act].split()[0])
                dacty=int(littlemat[act].split()[1])
                dactz=int(littlemat[act].split()[2])
                colidx = (topatom-1)*3 + (whichdir-1)
                #so 2  3  on first line means atom 2's z direction
                #then with atom 1x, 1y, 1z; 2x, 2y, 2z, etc.
                myhess[colidx][act*3+0]=dactx
                myhess[colidx][act*3+1]=dacty
                myhess[colidx][act*3+2]=dactz
                act = act + 1
            mycount = mycount + datoms + 1
        print(myhess)
        print numpy.linalg.eig(myhess)[0]
        return myhess
