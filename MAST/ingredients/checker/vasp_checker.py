from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio.vasp_output import Vasprun

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import dirutil
from MAST.utility.mastfile import MASTFile
from MAST.utility import MASTError
from MAST.utility.metadata import Metadata
import os
import shutil

def get_structure_from_file(filepath):
    """Get structure from file."""
    return Poscar.from_file(filepath, False).structure

def get_structure_from_directory(dirname):
    """Get structure from directory. Preferentially gets CONTCAR first."""
    cpath = os.path.join(dirname, "CONTCAR")
    ppath = os.path.join(dirname, "POSCAR")
    if os.path.isfile(cpath):
        return Poscar.from_file(cpath, False).structure
    elif os.path.isfile(ppath):
        return Poscar.from_file(ppath, False).structure
    else:
        raise MASTError("vasp_checker, get_structure_from_directory", "No valid structure in %s" % dirname)

def forward_parent_dynmat(parentpath, childpath, newname1="DYNMAT", newname2="XDATCAR"):
    """Forward the DYNMAT and the XDATCAR."""
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "DYNMAT"),os.path.join(childpath, newname1))
    shutil.copy(os.path.join(parentpath, "XDATCAR"),os.path.join(childpath, newname2))
    dirutil.unlock_directory(childpath)
    return


def forward_parent_structure(parentpath, childpath, newname="POSCAR"):
    """Copy CONTCAR to new POSCAR"""
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "CONTCAR"),os.path.join(childpath, newname))
    dirutil.unlock_directory(childpath)
    return

def forward_parent_initial_structure(parentpath, childpath, newname="POSCAR"):
    """Copy POSCAR to new POSCAR (for use in phonons)"""
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "POSCAR"),os.path.join(childpath, newname))
    dirutil.unlock_directory(childpath)
    return
def forward_parent_energy(parentpath, childpath, newname="OSZICAR"):
    """Copy OSZICAR"""
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
    """Check if single VASP non-NEB calculation is complete."""
    try:
        myoutcar = Outcar(os.path.join(dirname, "OUTCAR"))
    except (IOError):
        return False

    #hw 04/19/13
    try:
        if myoutcar.run_stats['User time (sec)'] > 0:
            return True
        else:
            return False
    except KeyError:
        return False

def is_ready_to_run(dirname):
    """Check if vasp is ready to run."""
    notready=0
    if not(os.path.isfile(dirname + "/KPOINTS")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/POTCAR")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/INCAR")):
        notready = notready + 1
    else: #INCAR exists. Now check POSCARs.
        myincar = MASTFile(dirname + "/INCAR")
        if myincar.get_line_match("IMAGES") == None: 
            if not(os.path.isfile(dirname + "/POSCAR")):
                notready = notready + 1
        else:
            subdirs = dirutil.walkdirs(dirname)
            subdirs.sort()
            if len(subdirs) == 0: #This is an NEB without folders yet
                notready = notready + 1
            else:
                for subdir in subdirs:
                    if not (os.path.isfile(subdir + "/POSCAR")):
                        notready = notready + 1
                if not os.path.isfile(subdirs[0] + "/OSZICAR"):
                    notready = notready + 1
                if not os.path.isfile(subdirs[-1] + "/OSZICAR"):
                    notready = notready + 1
    if not(os.path.isfile(dirname + "/submit.sh")):
        notready = notready + 1
    if notready > 0:
        return False
    else:
        return True

def _vasp_poscar_setup(keywords):
    name = keywords['name']
    pospath = os.path.join(name, "POSCAR")
    if os.path.isfile(pospath):
        my_poscar = Poscar.from_file(pospath) 
        #parent should have given a structure
    else: #this is an originating run; mast should give it a structure
        my_poscar = Poscar(keywords['structure'])
        dirutil.lock_directory(name)
        my_poscar.write_file(pospath)
        dirutil.unlock_directory(name)
    return my_poscar

def _vasp_kpoints_setup(keywords):
    """Parse mast_kpoints string, which should take the format:
        number, number, number designation
        examples: "3x3x3 M", "1x1x1 G". If no designation is given,
        Monkhorst-Pack is assumed.
    """
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
    name = keywords['name']

    if 'mast_xc' in keywords['program_keys'].keys():
        myxc = keywords['program_keys']['mast_xc'].upper() #Uppercase
    else:
        raise MASTError("vasp_checker, _vasp_potcar_setup","Exchange correlation functional needs to be specified in ingredients keyword mast_xc")

    if ('mast_pp_setup' in keywords['program_keys'].keys()):
        sites = my_poscar.site_symbols
        setup = keywords['program_keys']['mast_pp_setup']
        pp_sites = list()
        for site in sites:
            try:
                pp_sites.append(setup[site])
            except KeyError:
                pp_sites.append(site)
        my_potcar = Potcar(symbols=pp_sites, functional=myxc, sym_potcar_map=None)
    else:
        my_potcar = Potcar(symbols=my_poscar.site_symbols, functional=myxc, sym_potcar_map=None)

    dirutil.lock_directory(name)
    my_potcar.write_file(name + "/POTCAR")
    dirutil.unlock_directory(name)
    return my_potcar

def _vasp_incar_get_non_mast_keywords(program_keys_dict):
    """Get the non-VASP keywords and make a dictionary."""
    incar_dict=dict()
    allowedpath = os.path.join(dirutil.get_mast_install_path(), 'MAST',
                    'ingredients','programkeys','vasp_allowed_keywords.py')
    allowed_list = _vasp_incar_get_allowed_keywords(allowedpath)
    for key, value in program_keys_dict.iteritems():
        if not key[0:5] == "mast_":
            keytry = key.upper()
            if not (keytry in allowed_list):
                print "Ignoring program key %s for INCAR. To allow this keyword, add it to %s" % (keytry, allowedpath)
            else:
                if type(value)==str and value.isalpha():
                    incar_dict[keytry]=value.capitalize() #First letter cap
                else:
                    incar_dict[keytry]=value
    return incar_dict

def _vasp_is_neb(keywords):
    """Check if ingredient type is an actual NEB run.
        Returns:
            True if NEB, False otherwise.
    """
    metapath = "%s/metadata.txt" % keywords['name']
    if not os.path.isfile(metapath): #we are in some sort of sub-ingredient folder masquerading as a separate ingredient
        return False
    mymeta = Metadata(metafile="%s/metadata.txt" % keywords['name'])
    [ingline,ingval]=mymeta.search_data("ingredient type")
    if ingval == "NEB" or ingval == "NEBLowMesh":
        return True
    else:
        return False

def _vasp_incar_get_allowed_keywords(allowedpath):
    """Get allowed vasp keywords.
        Args:
            allowedpath <str>: file path for allowed vasp keywords
    """
    allowed = MASTFile(allowedpath)
    allowed_list=list()
    for entry in allowed.data:
        allowed_list.append(entry.strip())
    return allowed_list

def _vasp_incar_setup(keywords, my_potcar, my_poscar):
    """Set up the INCAR, including MAGMOM string, ENCUT, and NELECT."""
    name=keywords['name']
    myd = dict()
    myd = _vasp_incar_get_non_mast_keywords(keywords['program_keys'])
    if not _vasp_is_neb(keywords):
        try:
            myd.pop("IMAGES")
        except KeyError:
            pass
    if 'mast_multiplyencut' in keywords['program_keys'].keys():
        mymult = float(keywords['program_keys']['mast_multiplyencut'])
    else:
        mymult = 1.5
    if 'ENCUT' in myd.keys():
        pass
    else:
        myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(my_potcar)*mymult
    if 'mast_setmagmom' in keywords['program_keys'].keys():
        magstr = str(keywords['program_keys']['mast_setmagmom'])
        magmomstr=""
        maglist = magstr.split()
        numatoms = sum(my_poscar.natoms)
        if len(maglist) < numatoms:
            magct=0
            while magct < len(maglist):
                magmomstr = magmomstr + str(my_poscar.natoms[magct]) + "*" + maglist[magct] + " " 
                magct = magct + 1
        else:
            magmomstr = magstr
        myd['MAGMOM']=magmomstr
    if 'mast_charge' in keywords['program_keys'].keys():
        myelectrons = vasp_extensions.get_total_electrons(my_poscar, my_potcar)
        newelectrons=0.0
        try:
            adjustment = float(keywords['program_keys']['mast_charge'])
        except (ValueError, TypeError):
            raise MASTError("vasp_checker, vasp_incar_setup","Could not parse adjustment")
        #newelectrons = myelectrons + adjustment
        newelectrons = myelectrons - adjustment
        myd['NELECT']=str(newelectrons)
    my_incar = Incar(myd)
    dirutil.lock_directory(name)
    my_incar.write_file(name + "/INCAR")
    dirutil.unlock_directory(name)
    return my_incar

def set_up_program_input(keywords):
    myposcar = _vasp_poscar_setup(keywords)
    mykpoints = _vasp_kpoints_setup(keywords)
    mypotcar = _vasp_potcar_setup(keywords, myposcar)
    myincar = _vasp_incar_setup(keywords, mypotcar, myposcar)
    return

def get_path_to_write_neb_parent_energy(myname, myimages, parent):
    if parent == 1:
        return os.path.join(myname, "00", "OSZICAR")
    elif parent == 2:
        return os.path.join(myname, str(int(myimages)+1).zfill(2),"OSZICAR")
    elif len(parent) > 1:
        return os.path.join(myname, parent, "OSZICAR")
    else:
        raise MASTError("vasp_checker, get_path_to_write_neb_parent_energy","Parent not specified correctly.")

def set_up_neb_folders(myname, image_structures):
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
    set_up_neb_folders(keywords['name'], image_structures)
    mykpoints = _vasp_kpoints_setup(keywords)
    mypotcar = _vasp_potcar_setup(keywords, Poscar(image_structures[0]))
    myincar = _vasp_incar_setup(keywords, mypotcar, Poscar(image_structures[0]))
    return
def forward_extra_restart_files(parentpath, childpath):
    """Forward extra restart files: softlink to WAVECAR and CHGCAR."""
    dirutil.lock_directory(childpath)
    import subprocess
    os.chdir(childpath)
    print parentpath
    print childpath
    mylink=subprocess.Popen("ln -s %s/WAVECAR WAVECAR" % parentpath, shell=True)
    mylink.wait()
    mylink2=subprocess.Popen("ln -s %s/CHGCAR CHGCAR" % parentpath, shell=True)
    mylink2.wait()
    os.chdir(parentpath)
    dirutil.unlock_directory(childpath)
def add_selective_dynamics_to_structure(keywords, sdarray):
    name = keywords['name']
    pname = os.path.join(name,"POSCAR")
    phposcar = Poscar.from_file(pname)
    phposcar.selective_dynamics = sdarray
    dirutil.lock_directory(name)
    os.rename(pname, pname + "_no_sd")
    phposcar.write_file(pname)
    dirutil.unlock_directory(name)
    return

def get_vasp_energy(abspath):
    return Vasprun('%s/vasprun.xml' % abspath).ionic_steps[-1]["electronic_steps"][-1]["e_0_energy"]

