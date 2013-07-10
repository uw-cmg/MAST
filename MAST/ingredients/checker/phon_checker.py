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

def get_structure_from_file(filepath):
    """Get structure from file."""
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    return Poscar.from_file(filepath, False).structure

def get_structure_from_directory(dirname):
    """Get structure from directory. Preferentially gets CONTCAR first."""
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    cpath = os.path.join(dirname, "CONTCAR")
    ppath = os.path.join(dirname, "POSCAR")
    if os.path.isfile(cpath):
        return Poscar.from_file(cpath, False).structure
    elif os.path.isfile(ppath):
        return Poscar.from_file(ppath, False).structure
    else:
        raise MASTError("vasp_checker, get_structure_from_directory", "No valid structure in %s" % dirname)

def forward_parent_structure(parentpath, childpath, newname="POSCAR"):
    """Copy CONTCAR to new POSCAR"""
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "CONTCAR"),os.path.join(childpath, newname))
    dirutil.unlock_directory(childpath)
    return

def forward_parent_energy(parentpath, childpath, newname="OSZICAR"):
    """Copy OSZICAR"""
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "OSZICAR"),os.path.join(childpath, newname))
    dirutil.unlock_directory(childpath)
    return

def images_complete(dirname, numim):
    """Check if all images in a VASP NEB calculation are complete.
        dirname = directory housing /00.../0N+1 files; 
                  only checks directories /01.../0N where N is # images
        numim = number of images
    """
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    imct=1
    numim = int(numim)
    while imct <= numim:
        num_str = str(imct).zfill(2)
        impath = os.path.join(dirname, num_str)
        try:
            myoutcar = Outcar(os.path.join(impath, "OUTCAR"))
        except (IOError):
            return False
        if not 'User time (sec)' in myoutcar.run_stats.keys():
            return False
        if myoutcar.run_stats['User time (sec)'] > 0:
            #print "image",imct,"complete"
            pass
        else:
            return False
        imct = imct + 1
    return True


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

def _vasp_kpoints_setup(keywords):
    """Parse mast_kpoints string, which should take the format:
        number, number, number designation
        examples: "3x3x3 M", "1x1x1 G". If no designation is given,
        Monkhorst-Pack is assumed.
    """
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    name = keywords['name']
    if 'mast_kpoints' in keywords['program_keys'].keys():
        kptlist = keywords['program_keys']['mast_kpoints']
    else:
        raise MASTError("vasp_checker, _vasp_kpoint_setup","k-point instructions need to be set in ingredients keyword mast_kpoints")
    if len(kptlist) == 3:
        desig = "M"
    else:
        desig = kptlist[3].upper()
    if desig == "M":
        my_kpoints = Kpoints.monkhorst_automatic(kpts=(int(kptlist[0]),int(kptlist[1]),int(kptlist[2])),shift=(0,0,0))
    elif desig == "G":
        my_kpoints = Kpoints.gamma_automatic(kpts=(int(kptlist[0]),int(kptlist[1]),int(kptlist[2])),shift=(0,0,0))
    else:
        raise MASTError("vasp_checker, _vasp_kpoint_setup","kpoint designation " + desig + " not recognized.")
    dirutil.lock_directory(name)
    my_kpoints.write_file(name + "/KPOINTS")
    dirutil.unlock_directory(name)
    return my_kpoints

def _vasp_potcar_setup(keywords, my_poscar):
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    name=keywords['name']
    if 'mast_xc' in keywords['program_keys'].keys():
        myxc = keywords['program_keys']['mast_xc'].upper() #Uppercase
    else:
        raise MASTError("vasp_checker, _vasp_potcar_setup","Exchange correlation functional needs to be specified in ingredients keyword mast_xc")
    my_potcar = Potcar(symbols=my_poscar.site_symbols, functional=myxc, sym_potcar_map=None)
    dirutil.lock_directory(name)
    my_potcar.write_file(name + "/POTCAR")
    dirutil.unlock_directory(name)
    return my_potcar

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
    _nosd_my_dynmat(keywords)
    name=keywords['name']
    if not os.path.isfile(name + "/DYNMAT_mod"):
        raise MASTError("checker/phon_checker", "No modified DYNMAT found in %s." % name)
    myforces=MASTFile(name + "/DYNMAT_mod")
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
    if not os.path.isfile(name + "/DYNMAT"):
        raise MASTError("checker/phon_checker", "No DYNMAT found in %s." % name)
    myforces=MASTFile(name + "/DYNMAT")
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
    myforces.to_file(name + "/DYNMAT_mod")



def set_up_program_input(keywords):
    _phon_poscar_setup(keywords)
    _phon_inphon_setup(keywords)
    _phon_forces_setup(keywords)
    return

def get_path_to_write_neb_parent_energy(myname, myimages, parent):
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    if parent == 1:
        return os.path.join(myname, "00", "OSZICAR")
    elif parent == 2:
        return os.path.join(myname, str(int(myimages)+1).zfill(2),"OSZICAR")
    elif len(parent) > 1:
        return os.path.join(myname, parent, "OSZICAR")
    else:
        raise MASTError("vasp_checker, get_path_to_write_neb_parent_energy","Parent not specified correctly.")

def set_up_neb_folders(myname, image_structures):
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    imct=0
    while imct < len(image_structures):
        imposcar = Poscar(image_structures[imct])
        num_str = str(imct).zfill(2)
        impath = os.path.join(myname, num_str)
        impospath = os.path.join(myname, "POSCAR_" + num_str)
        dirutil.lock_directory(myname)
        imposcar.write_file(impospath)
        dirutil.unlock_directory(myname)
        try:
            os.makedirs(impath)
        except OSError:
            print "Directory at", impath, "already exists."
            return None
        dirutil.lock_directory(impath)
        imposcar.write_file(os.path.join(impath, "POSCAR"))
        dirutil.unlock_directory(impath)
        imct = imct + 1
    return
    

def set_up_program_input_neb(keywords, image_structures):
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    set_up_neb_folders(keywords['name'], image_structures)
    mykpoints = _vasp_kpoints_setup(keywords)
    mypotcar = _vasp_potcar_setup(keywords, Poscar(image_structures[0]))
    myincar = _vasp_incar_setup(keywords, mypotcar, Poscar(image_structures[0]))
    return

def add_selective_dynamics_to_structure(keywords, sdarray):
    raise MASTError("checker/phon_checker","Not implemented for PHON!")
    name = keywords['name']
    pname = os.path.join(name,"POSCAR")
    phposcar = Poscar.from_file(pname)
    phposcar.selective_dynamics = sdarray
    dirutil.lock_directory(name)
    os.rename(pname, pname + "_no_sd")
    phposcar.write_file(pname)
    dirutil.unlock_directory(name)
    return
