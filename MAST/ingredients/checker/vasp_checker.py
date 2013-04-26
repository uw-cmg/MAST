from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from MAST.utility import dirutil
import os
import shutil

def get_structure_from_file(filepath):
    """Get structure from file."""
    return Poscar.from_file(filepath, False).structure

def forward_parent_structure(parentpath, childpath, newname="POSCAR"):
    """Copy CONTCAR to new POSCAR"""
    dirutil.lock_directory(childpath)
    shutil.copy(os.path.join(parentpath, "CONTCAR"),os.path.join(childpath, newname))
    dirutil.unlock_directory(childpath)
    return

def images_complete(dirname, numim):
    """Check if all images in a VASP NEB calculation are complete.
        dirname = directory housing /00.../0N+1 files; 
                  only checks directories /01.../0N where N is # images
        numim = number of images
    """
    imct=1
    while imct <= numim:
        num_str = str(imct).zfill(2)
        impath = os.path.join(dirname, num_str)
        try:
            myoutcar = Outcar(os.path.join(impath, "OUTCAR"))
        except (IOError):
            return False
        if myoutcar.run_stats['User time (sec)'] > 0:
            pass
        else:
            return False
        imct = imct + 1
    return True


def is_complete(dirname):
    """Check if all images in a VASP NEB calculation are complete.
        dirname = directory housing /00.../0N+1 files; 
                  only checks directories /01.../0N where N is # images
        numim = number of images
    """
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
    if not(os.path.isfile(dirname + "/POSCAR")):
        notready = notready + 1
    if not(os.path.isfile(dirname + "/submit.sh")):
        notready = notready + 1
    if notready > 0:
        return False
    else:
        return True
