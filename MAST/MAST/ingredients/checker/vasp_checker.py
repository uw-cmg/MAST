from pymatgen.io.vaspio import Poscar
import os


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

