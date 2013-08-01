from pymatgen.io.vaspio import *
import numpy as np
from MAST.utility.dirutil import *
from MAST.utility import MASTError
from MAST.utility import MASTFile
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

def read_my_dynmat(mydir, fname="DYNMAT"):
    """Read a DYNMAT file.
        Returns:
            dyndict <dict>: dictionary structured like this:
                ['numspec'] = <int> number of species
                ['numatoms'] = <int> number of atoms
                ['numdisp'] = <int> number of displacements
                ['massline'] = <str> masses line
                ['atoms'][atom <int>][disp <int>]['displine'] = 
                    displacement
                    vector (part of first line in dynmat block, 
                    e.g. "0.01 0 0")
                ['atoms'][atom <int>][disp <int>]['dynmat'] = 
                        <list> 
                        list of dynmat lines for this atom 
                        and this displacement
    """
    mydyn = MASTFile(os.path.join(mydir, fname))
    mydata = list(mydyn.data) #whole new copy
    dyndict=dict()
    firstline = mydata.pop(0) #pop first value
    firstspl = firstline.strip().split()
    dyndict['numspec']=int(firstspl[0])
    numatoms = int(firstspl[1])
    dyndict['numatoms']=numatoms
    dyndict['numdisp']=int(firstspl[2])
    dyndict['massline']=mydata.pop(0)
    dyndict['atoms']=dict()
    atom=0
    disp=0
    thirdline=""
    while len(mydata) > 0:
        thirdline = mydata.pop(0)
        thirdspl = thirdline.strip().split()
        atom=int(thirdspl[0])
        disp=int(thirdspl[1])
        displine=' '.join(thirdspl[2:])
        if not atom in dyndict['atoms'].keys():
            dyndict['atoms'][atom]=dict()
        dyndict['atoms'][atom][disp]=dict()
        dyndict['atoms'][atom][disp]['displine'] = displine
        dyndict['atoms'][atom][disp]['dynmat']=list()
        for act in range(0, numatoms):
            dyndict['atoms'][atom][disp]['dynmat'].append(mydata.pop(0))
    return dyndict

def write_my_dynmat(mydir, dyndict, fname="DYNMAT"):
    """Write a DYNMAT file.
        Args:
            mydir <str>: Directory in which to write
            dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
            fname <str>: filename (default DYNMAT)
    """
    dynwrite=MASTFile()
    dynwrite.data=list()
    firstline=str(dyndict['numspec']) + " " + str(dyndict['numatoms']) + " " + str(dyndict['numdisp']) + "\n"
    dynwrite.data.append(firstline)
    dynwrite.data.append(dyndict['massline'])
    atomlist=dyndict['atoms'].keys()
    atomlist.sort()
    for atom in atomlist:
        displist = dyndict['atoms'][atom].keys()
        displist.sort()
        for disp in displist:
            thirdline = str(atom) + " " + str(disp) + " " + dyndict['atoms'][atom][disp]['displine'] + "\n"
            dynwrite.data.append(thirdline)
            for line in dyndict['atoms'][atom][disp]['dynmat']:
                dynwrite.data.append(line)
    dynwrite.to_file(os.path.join(mydir, fname))
    

def write_my_dynmat_without_disp_or_mass(mydir, dyndict, fname="DYNMAT"):
    """Write a DYNMAT file without the displacement indicators 1, 2, 3 
        and without the masses line, and with first line having only
        the total number of displacements, for PHON.
        Args:
            mydir <str>: Directory in which to write
            dyndict <dict>: Dictionary of dynmat (see read_my_dynmat)
            fname <str>: filename (default DYNMAT)
    """
    dynwrite=MASTFile()
    dynwrite.data=list()
    firstline=str(dyndict['numdisp']) + "\n"
    dynwrite.data.append(firstline)
    atomlist=dyndict['atoms'].keys()
    atomlist.sort()
    for atom in atomlist:
        displist = dyndict['atoms'][atom].keys()
        displist.sort()
        for disp in displist:
            thirdline = str(atom) + " " + dyndict['atoms'][atom][disp]['displine'] + "\n"
            dynwrite.data.append(thirdline)
            for line in dyndict['atoms'][atom][disp]['dynmat']:
                dynwrite.data.append(line)
    dynwrite.to_file(os.path.join(mydir, fname))

def read_my_xdatcar(mydir, fname="XDATCAR"):
    """Read an XDATCAR file.
        Returns: 
            xdatdict <dict>: Dictionary of configurations
                xdatdict['descline'] = <str> description line
                xdatdict['specline'] = <str> species line
                xdatdict['numline'] = <str> numbers lines
                xdatdict['type'] = <str> Direct or not
                xdatdict['numatoms'] = <int> number of atoms
                xdatdict['configs'][config <int>] = <list>
                        list of lines in this configuration
    """
    myxdat = MASTFile(os.path.join(mydir, fname))
    mydata = list(myxdat.data) #whole new copy
    xdatdict=dict()
    xdatdict['descline'] = mydata.pop(0) #pop first value
    xdatdict['specline'] = mydata.pop(0)
    numline = mydata.pop(0)
    xdatdict['numline'] = numline
    numatoms = sum(map(int, numline.strip().split()))
    xdatdict['numatoms'] = numatoms
    xdatdict['type'] = mydata.pop(0)
    xdatdict['configs']=dict()
    kfgct=1
    while len(mydata) > 0:
        mydata.pop(0) # Konfig line, or a blank line
        xdatdict['configs'][kfgct] = list()
        for act in range(0,numatoms):
            xdatdict['configs'][kfgct].append(mydata.pop(0))
        kfgct=kfgct+1
    return xdatdict

def write_my_xdatcar(mydir, xdatdict, fname="XDATCAR"):
    """Write an XDATCAR file.
        Args:
            mydir <str>: Directory in which to write
            xdatdict <dict>: Dictionary of XDATCAR (see read_my_xdatcar)
            fname <str>: filename (default DYNMAT)
    """
    xdatwrite=MASTFile()
    xdatwrite.data=list()
    xdatwrite.data.append(xdatdict['descline'])
    xdatwrite.data.append(xdatdict['specline'])
    xdatwrite.data.append(xdatdict['numline'])
    xdatwrite.data.append(xdatdict['type'])
    configlist = xdatdict['configs'].keys()
    configlist.sort()
    for cfg in configlist:
        xdatwrite.data.append("Konfig=%1i\n" % cfg)
        xdatwrite.data.extend(xdatdict['configs'][cfg])
    xdatwrite.to_file(os.path.join(mydir, fname))

def combine_dynmats(mydir):
    """Combine DYNMATs into one file.
        Args:
            mydir <str>: top directory for DYNMAT files
    """
    dynmatlist = walkfiles(mydir, 2, 5, "*DYNMAT*") #start one level below
    if len(dynmatlist) == 0:
        raise MASTError("pmgextend combine_dynmats", "No DYNMATs found under " + mydir)
    totnumdisp=0
    largedyn=dict()
    largedyn['atoms'] = dict()
    for onedynmat in dynmatlist:
        dyndir = os.path.dirname(onedynmat)
        onedyn = read_my_dynmat(dyndir)
        totnumdisp = totnumdisp + onedyn['numdisp']
        for atom in onedyn['atoms'].keys():
            if not atom in largedyn['atoms'].keys():
                largedyn['atoms'][atom]=dict()
                mydisp=1 #start at 1
            for disp in onedyn['atoms'][atom].keys():
                if disp in largedyn['atoms'][atom].keys():
                    mydisp=mydisp + 1 #increment
                largedyn['atoms'][atom][mydisp]=dict()
                largedyn['atoms'][atom][mydisp]['displine'] = str(onedyn['atoms'][atom][disp]['displine'])
                largedyn['atoms'][atom][mydisp]['dynmat']=list(onedyn['atoms'][atom][disp]['dynmat'])
    largedyn['numspec'] = onedyn['numspec'] #should all be the same
    largedyn['numatoms'] = onedyn['numatoms']
    largedyn['massline'] = onedyn['massline']
    largedyn['numdisp'] = totnumdisp
    write_my_dynmat(mydir, largedyn, "DYNMAT_combined")
    
def combine_displacements(mydir):
    """Combine displacements (here XDATCARs) into one file.
        Args:
            mydir <str>: top directory for DYNMAT files
    """
    largexdat=dict()
    xdatlist = walkfiles(mydir, 2, 5, "*XDATCAR*") #start one level below
    if len(xdatlist) == 0:
        raise MASTError("pmgextend combine_displacements", "No XDATCARs found under " + mydir)
    configs=list()
    kfgct=1 # skip config 1 until the end
    largexdat['configs']=dict()
    xdatlist.sort() #get them in order
    for onexdatmat in xdatlist:
        xdatdir = os.path.dirname(onexdatmat)
        onexdat = read_my_xdatcar(xdatdir)
        for kfg in onexdat['configs'].keys():
            kfgct = kfgct + 1 #start at 2
            if kfg == 1: #skip config 1 until the end
                pass
            else:
                largexdat['configs'][kfgct] = onexdat['configs'][kfg]
    largexdat['configs'][1] = onexdat['configs'][1] #get one of the first configs (the first should be the same for all the XDATCARs)
    largexdat['descline'] = onexdat['descline'] #should all be the same
    largexdat['specline'] = onexdat['specline']
    largexdat['numline'] = onexdat['numline']
    largexdat['numatoms'] = onexdat['numatoms']
    largexdat['type'] = onexdat['type']
    write_my_xdatcar(mydir, largexdat, "XDATCAR_combined")
def make_hessian(myposcar, mydir):
    """Combine DYNMATs into one hessian and solve for frequencies.
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

