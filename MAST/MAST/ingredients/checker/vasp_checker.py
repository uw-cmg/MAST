from pymatgen.io.vaspio import Poscar
import os


def get_structure_from_parent(parentpath):
    """Get structure from VASP output."""
    structpath = os.path.join(parentpath, "CONTCAR")
    parent_structure = Poscar.from_file(structpath, False).structure
    return parent_structure

