from pymatgen.io.vaspio import *

def get_max_enmax_from_Potcar(mypotcar):
    """Get maximum enmax value (float) from Potcar (combined list)"""
    enmax_list=list()
    potcarct=0
    onepotcar=None
    while potcarct < len(mypotcar):
        onepotcar = mypotcar[potcarct] #A PotcarSingle object
        enmax_list.append(onepotcar.enmax)
        potcarct = potcarct + 1
    return max(enmax_list)

def make_one_unfrozen_atom_poscar(myposcar, natom)
    mysd=np.zeros([sum(myposcar.natoms),3],bool)
    mysd[natom][0]=True
    mysd[natom][1]=True
    mysd[natom][2]=True
    myposcar.selective_dynamics = mysd
    return myposcar

