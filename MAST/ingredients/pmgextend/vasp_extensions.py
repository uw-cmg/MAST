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

    """Recombine DYNMAT files from one unfrozen atom in folder."""
    natoms = sum(myposcar.natoms)
    myhess=np.zeros([natoms*3, natoms*3])
    #arrange as x1 y1 z1 x2 y2 z2 x3 y3 z3...
    dirlist = os.listdir(mydir)
    dynlist=[]
#    for entry in dirlist:

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
    print "DYNMATLIST:"
    print dynmatlist
    datoms=0
    dmats=0
    mct=0
    for onedynmat in dynmatlist:
        dynlines=[]
        opendyn = open(onedynmat,'rb')
        dynlines=opendyn.readlines()
        opendyn.close()
        datoms = int(dynlines[0].split()[1])
        dmats = int(dynlines[0].split()[2])
        mycount=2 #starting line
        while mycount < len(dynlines)-datoms:
            littlemat=[]
            topatom = int(dynlines[mycount].split()[0])
            whichdir = int(dynlines[mycount].split()[1])
            littlemat = dynlines[mycount+1:mycount+datoms+1]
            print littlemat
            act = 0
            while act < datoms:
                dactx=float(littlemat[act].split()[0])
                dacty=float(littlemat[act].split()[1])
                dactz=float(littlemat[act].split()[2])
                colidx = (topatom-1)*3 + (whichdir-1)
                #so 2  3  on first line means atom 2's z direction
                #then with atom 1x, 1y, 1z; 2x, 2y, 2z, etc.
                myhess[colidx][act*3+0]=dactx
                myhess[colidx][act*3+1]=dacty
                myhess[colidx][act*3+2]=dactz
                act = act + 1
            mycount = mycount + datoms + 1
            print mycount
    print "UNALTERED HESSIAN:"
    print(myhess)
    #create mass matrix
    masses=dynlines[1].split()
    print "MASSES:", masses
    massarr=np.zeros([datoms*3,1])
    act=0
    print myposcar.natoms
    nspec=len(myposcar.natoms)
    totatoms=0
    while act < datoms:
        mymass=0
        nct=0
        totatoms=0
        while (mymass==0) and nct < nspec:
            totatoms = totatoms + myposcar.natoms[nct]
            if act < totatoms:
                mymass = float(masses[nct])
            nct = nct + 1
        print mymass
        massarr[act*3+0][0]=mymass
        massarr[act*3+1][0]=mymass
        massarr[act*3+2][0]=mymass
        act = act + 1
    massmat = massarr*np.transpose(massarr)
    print "MASS MAT:"
    print massmat
    print "STEP:"
    step = float(dynlines[2].split()[2])
    print step
    print "HESSIAN * -1 / step / sqrt(mass1*mass2)" 
    normhess=np.zeros([natoms*3,natoms*3])
    cidx=0
    while cidx < natoms*3:
        ridx=0
        while ridx < natoms*3:
            normhess[ridx][cidx]=-1*myhess[ridx][cidx]/step/np.sqrt(massmat[ridx][cidx])
            ridx = ridx + 1
        cidx = cidx + 1
    print normhess
    print "EIGENVALUES:"
    myeig = np.linalg.eig(normhess)[0]
    print myeig
    print "SQRT of EIGENVALUES in sqrt(eV/AMU)/Angstrom/2pi:"
    myfreq = np.sqrt(myeig)
    print myfreq
    print "SQRT OF EIGENVALUES in THz:"
    myfreqThz = myfreq*15.633302
    print myfreqThz
    myfreqThzsorted = myfreqThz
    myfreqThzsorted.sort()
    print myfreqThzsorted
    return myfreqThzsorted

def get_total_electrons(myposcar, mypotcar):
    """Get the total number of considered electrons in the system."""
    atomlist = myposcar.natoms
    zvallist = get_zval_list(mypotcar)
    totzval = 0.0
    atomct = 0
    if not (len(zvallist) == len(atomlist)):
        raise MASTError("pmgextend, get_total_electrons",
            "Number of species and number of POTCARs do not match.")
    while atomct < len(atomlist):
        totzval = totzval + (atomlist[atomct] * zvallist[atomct])
        atomct = atomct + 1
    return totzval

def get_zval_list(mypotcar):
    """Get zvals from POTCAR"""
    zval_list=list()
    potcarct=0
    onepotcar=None
    while potcarct < len(mypotcar):
        onepotcar = mypotcar[potcarct] #A PotcarSingle object
        zval_list.append(onepotcar.zval)
        potcarct = potcarct + 1
    return zval_list

