from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar
from pymatgen.io.vaspio import Kpoints
from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
import os
import shutil
import pymatgen
import numpy as np


def is_complete(dirname):
    """Check if PHON thermo run is complete."""
    if os.path.isfile(os.path.join(dirname, "THERMO")):
        return True
    else:
        return False

def is_ready_to_run(dirname):
    """Check if PHON is ready to run."""
    notready=0
    if not(os.path.isfile(dirname + "/FORCES")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/POSCAR")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/INPHON")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/submit.sh")):
        notready = notready + 1
    if notready > 0:
        return False
    else:
        return True

def _phon_poscar_setup(keywords):
    """Set up a PHON POSCAR file. Strip out the "elements" line (that is,
        use VASP version 4 format. Also strip out anything beneath the atoms
        line.
    """
    name = keywords['name']
    pospath = os.path.join(name, "POSCAR")
    if os.path.isfile(pospath + "_prePHON"): #Already done. Return.
        return
    if os.path.isfile(pospath):
        my_poscar = Poscar.from_file(pospath) 
        #parent should have given a structure
    else: #this is an originating run; mast should give it a structure
        my_poscar = Poscar(keywords['structure'])
    my_poscar.selective_dynamics=None #unset SD if it is set
    my_poscar.velocities=None #unset velocities
    #write two copies
    dirutil.lock_directory(name)
    my_poscar.write_file(pospath)
    my_poscar.write_file(pospath + "_prePHON")
    dirutil.unlock_directory(name)
    #pick up a copy and strip out the elements line.
    mypfile = MASTFile(pospath)
    myline6=mypfile.get_line_number(6)
    if myline6.strip().split()[0].isalpha:
        mypfile.modify_file_by_line_number(6,"D")
    mypfile.to_file(pospath)
    return



def _phon_inphon_get_non_mast_keywords(program_keys_dict):
    """Sort out the non-PHON keywords and make a dictionary."""
    inphon_dict=dict()
    allowedpath = os.path.join(dirutil.get_mast_install_path(), 'MAST',
                    'ingredients','programkeys','phon_allowed_keywords.py')
    allowed_list = _phon_inphon_get_allowed_keywords(allowedpath)
    for key, value in program_keys_dict.iteritems():
        if not key[0:5] == "mast_":
            keytry = key.upper()
            if not (keytry in allowed_list):
                print "Ignoring program key %s for INPHON. To allow this keyword, add it to %s" % (keytry, allowedpath)
            else:
                if type(value)==str and value.isalpha():
                    inphon_dict[keytry]=value.capitalize() #First letter cap
                else:
                    inphon_dict[keytry]=value
    return inphon_dict

def _phon_inphon_get_allowed_keywords(allowedpath):
    """Get allowed PHON keywords.
        Args:
            allowedpath <str>: file path for allowed PHON keywords
    """
    allowed = MASTFile(allowedpath)
    allowed_list=list()
    for entry in allowed.data:
        allowed_list.append(entry.strip())
    return allowed_list

def _phon_inphon_setup(keywords):
    """Set up the INPHON file."""
    name=keywords['name']
    myd = dict()
    myd = _phon_inphon_get_non_mast_keywords(keywords['program_keys'])
    my_inphon = MASTFile()
    for key, value in myd.iteritems():
        my_inphon.data.append(str(key) + "=" + str(value).upper() + "\n")
    if not ("NTYPES" in myd.keys()) and not ("MASS" in myd.keys()):
        [nline,massline] = _phon_inphon_get_masses(keywords)
        my_inphon.data.append(nline + "\n")
        my_inphon.data.append(massline + "\n")
    my_inphon.to_file(name + "/INPHON")
    return 

def _phon_inphon_get_masses(keywords):
    """Get the ntypes and masses line for INPHON.
        Returns:
            nline, massline
            nline <str>: NYTPES = <number of atomic types in POSCAR>
            massline <str>: MASS = <mass1> <mass2> ...
    """
    name=keywords['name']
    if os.path.isfile(name + "/POSCAR_prePHON"):
        mypos = Poscar.from_file(name + "/POSCAR_prePHON")
    else: #not yet modified to strip out the species line.
        mypos = Poscar.from_file(name + "/POSCAR")
    sitesym = mypos.site_symbols
    nline="NTYPES=" + str(len(sitesym))
    massline="MASS="
    for sym in sitesym:
        el = pymatgen.core.periodic_table.Element(sym)
        massline = massline + str(el.atomic_mass) + " "
    return [nline, massline]

def _phon_forces_setup(keywords):
    """Set up the FORCES file. This is like the DYNMAT but with the mass
        line stripped out and no direction indicators. Also, a block must
        be present for every atom, with a displacement, even if all entries
        are zero (e.g. fake block for selective dynamics)
    """
    _replace_my_displacements(keywords)
    _nosd_my_dynmat(keywords)
    name=keywords['name']
    if not os.path.isfile(name + "/DYNMAT_mod_2"):
        raise MASTError("checker/phon_checker", "No DYNMAT_mod_2 found in %s." % name)
    myforces=MASTFile(name + "/DYNMAT_mod_2")
    infosplit = myforces.get_line_number(1).strip().split()
    numdisp = int(infosplit[2])  #number of dynmat chunks
    numatoms = int(infosplit[1]) #number of lines in a dynmat chunk
    idxrg=list()
    for nct in range(0,numdisp):
        idxrg.append(1+1+1+nct*(numatoms + 1)) #get all the header lines
    for idx in idxrg: #strip out the 'direction' indicator 1, 2, 3
        mysplit = myforces.get_line_number(idx).strip().split()
        mynewlist = list()
        mynewlist.append(mysplit[0])
        mynewlist.extend(mysplit[2:])
        mynewline = ' '.join(mynewlist) + "\n"
        myforces.modify_file_by_line_number(idx, "R", mynewline)
    myforces.modify_file_by_line_number(1, "R", str(numdisp) + "\n") #modify info line
    myforces.modify_file_by_line_number(2, "D") #remove masses line
    myforces.to_file(name + "/FORCES")
    return

def _nosd_my_dynmat(keywords):
    """Creates fake blocks in DYNMAT for filling back in the atoms and 
        directions skipped through selective dynamics.
    """
    name=keywords['name']
    if not os.path.isfile(name + "/DYNMAT_mod_1"):
        raise MASTError("checker/phon_checker", "No DYNMAT_mod_1 found in %s." % name)
    myforces=MASTFile(name + "/DYNMAT_mod_1")
    infosplit = myforces.get_line_number(1).strip().split()
    numspec = int(infosplit[0])
    numdisp = int(infosplit[2])  #number of dynmat chunks
    numatoms = int(infosplit[1]) #number of lines in a dynmat chunk
    idxrg=list()
    atomsanddirs=dict()
    for nct in range(0,numdisp):
        idxrg.append(1+1+1+nct*(numatoms + 1)) #get all the header lines
    for idx in idxrg: #strip out the 'direction' indicator 1, 2, 3
        mysplit = myforces.get_line_number(idx).strip().split()
        myatom = int(mysplit[0])
        if not myatom in atomsanddirs.keys():
            atomsanddirs[myatom] = dict()
            atomsanddirs[myatom][1] = 'not found'
            atomsanddirs[myatom][2] = 'not found'
            atomsanddirs[myatom][3] = 'not found'
        mydirection = int(mysplit[1])
        atomsanddirs[myatom][mydirection] = 'found'
    lct=2 #start at line 2 ("masses" line)
    zeroline="0 0 0" + "\n"
    for act in range(1,numatoms + 1):
        if not act in atomsanddirs.keys():
            atomsanddirs[act] = dict()
            atomsanddirs[act][1] = 'not found'
            atomsanddirs[act][2] = 'not found'
            atomsanddirs[act][3] = 'not found'
        for didx in range(1,4):
            if atomsanddirs[act][didx] == 'not found':
                if didx == 1:
                    headerline=str(act) + " 1 0.00001 0 0" + "\n"
                elif didx == 2:
                    headerline=str(act) + " 2 0 0.00001 0" + "\n"
                elif didx == 3:
                    headerline=str(act) + " 3 0 0 0.00001" + "\n"
                myforces.modify_file_by_line_number(lct, "I", headerline)
                lct = lct + 1
                for zct in range(0,numatoms):
                    myforces.modify_file_by_line_number(lct, "I", zeroline)
                    lct = lct + 1
            else:
                lct = lct + 1 + numatoms
    numdisp = numatoms*3
    newheader = str(numspec) + " " + str(numatoms) + " " + str(numdisp) + "\n"
    myforces.modify_file_by_line_number(1, "R", newheader)
    myforces.to_file(name + "/DYNMAT_mod_2")

def _replace_my_displacements(keywords):
    """
        In VASP, 0.01 0 0 in DYNMAT from phonons is 1/Angstroms in the 
        x-direction (XDATCAR shows that it makes fractional coord 
        displacements for all 3 lattice vectors in a non-cubic system to get 
        this strictly-x-direction) 
        In PHON, 0.01 0 0 means 0.01 multiplied by lattice vector a.
        Back out the fractional displacements used by VASP from the XDATCAR, 
        match them up, and use them.
        Konfig =1 is the un-displaced cell.
        Now for NFREE=2,
        there are two Konfigs for each displacment; the first is positive
        POTIM in the x-direction (for example, POTIM = 0.01), then negative
        POTIM in the x-direction, then y, then z.
        So one unfrozen atom has seven Konfigs.
        DYNMAT, however, reports the AVERAGE force from each Konfig pair.
        So we only want Konfigs 2, 4, and 6, corresponding to POTIM 0 0, 
        0 POTIM 0, and 0 0 POTIM
    """
    name=keywords['name']
    if not os.path.isfile(name + "/XDATCAR"):
        raise MASTError("checker/phon_checker", "No XDATCAR found in %s." % name)
    myxdat=vasp_extensions.read_my_xdatcar(name)
    if not os.path.isfile(name + "/DYNMAT"):
        raise MASTError("checker/phon_checker", "No DYNMAT found in %s." % name)
    myforces=vasp_extensions.read_my_dynmat(name)
    #resume editing here
    infosplit = myforces.get_line_number(1).strip().split()
    numspec = int(infosplit[0])
    numdisp = int(infosplit[2])  #number of dynmat chunks
    numatoms = int(infosplit[1]) #number of lines in a dynmat chunk
    idxrg=list() #list of header line numbers in DYNMAT
    kfgrg=list() #list of Konfigs in XDATCAR
    atomsanddirs=dict()
    for nct in range(0,numdisp*2+1): #get all konfig lines
        kfgrg.append(4+1+nct*(numatoms + 1))
    kfgdict=dict()
    kfgnum=0
    for kidx in kfgrg:
        if kidx == 5: #first configuration
            kct = 0
        if np.mod(kfgnum,2) == 1: #skip odd configurations
            pass
        else:
            kct = np.divide(kfgnum,2)
            kfgdict[kct]=dict()
            for act in range(1,numatoms+1):
                klin=kidx+act
                mycoords=np.array(myxdat.get_line_number(klin).strip().split(),float)
                kfgdict[kct][act]=mycoords
        kfgnum = kfgnum + 1
    import pdb
    pdb.set_trace()
    for nct in range(0,numdisp):
        idxrg.append(1+1+1+nct*(numatoms + 1)) #get all the header lines
    dct=1
    for idx in idxrg:
        mysplit = myforces.get_line_number(idx).strip().split()
        myatom = mysplit[0]
        mydispdir = mysplit[1]
        newcoords = kfgdict[dct][int(myatom)] - kfgdict[0][int(myatom)]
        coordstr = str(newcoords[0]) + " " + str(newcoords[1]) + " " + str(newcoords[2])
        mynewline = myatom + " " + mydispdir + " " + coordstr + "\n"
        myforces.modify_file_by_line_number(idx, "R", mynewline)
        dct=dct+1
    myforces.to_file(name + "/DYNMAT_mod_1")

def set_up_program_input(keywords):
    _phon_poscar_setup(keywords)
    _phon_inphon_setup(keywords)
    _phon_forces_setup(keywords)
    return

