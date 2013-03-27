from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
import os
import shutil

def get_structure_from_parent(parentpath):
    """Get structure from VASP output."""
    structpath = os.path.join(parentpath, "CONTCAR")
    try:
        parent_structure = Poscar.from_file(structpath, False).structure
    except:
        print "Error getting structure from CONCAR."
        print "Make except statement more specific."
        return None
    return parent_structure

def forward_parent_structure(parentpath, childpath):
    """Copy CONTCAR to new POSCAR"""
    shutil.copy(os.path.join(parentpath, "CONTCAR"),os.path.join(childpath, "POSCAR"))
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


