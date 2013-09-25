import os
import subprocess
import time

from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import Metadata
allowed_keys = {
    'name' : (str, str(), 'Name of directory'),
    'program': (str, str(), 'Program, e.g. "vasp"'),
    'program_keys': (dict, dict(), 'Dictionary of program keywords'),
    'structure': (Structure, None, 'Pymatgen Structure object')
            }
class BaseChecker(MASTObj):
    """Base checker class. This class switches between
        program-specific functions.
        For example, phon-related functions or vasp-related
        functions get their own checker class.
    """
    
    def __init__(self, allowed_keys, **kwargs):
        allowed_keys_base = dict()
        allowed_keys_base.update(allowed_keys) 
        MASTObj.__init__(self, allowed_keys_base, **kwargs)

    def get_structure_from_directory(self, dirname):
        raise NotImplementedError

    def get_structure_from_file(self, filepath):
        raise NotImplementedError

    def forward_parent_structure(self, parentpath, childpath, newname="POSCAR"):
        raise NotImplementedError
    
    def forward_parent_initial_structure(self, parentpath, childpath, newname="POSCAR"):
        raise NotImplementedError

    def forward_parent_energy(self, parentpath, childpath, newname="OSZICAR"):
        raise NotImplementedError
    def forward_parent_dynmat(self, parentpath, childpath, newname="DYNMAT"):
        raise NotImplementedError
    def is_complete(self):
        raise NotImplementedError
   
    def is_ready_to_run(self):
        raise NotImplementedError
    
    def set_up_program_input(self):
        raise NotImplementedError

    def get_path_to_write_neb_parent_energy(self, parent):
        raise NotImplementedError

    def set_up_program_input_neb(self, image_structures):
        raise NotImplementedError

    def add_selective_dynamics_to_structure(self, sdarray):
        raise NotImplementedError

    def forward_extra_restart_files(self, parentpath, childpath):
        raise NotImplementedError
    
    def combine_dynmats(self):
        raise NotImplementedError
    
    def combine_displacements(self):
        raise NotImplementedError
    
    def get_e0_energy(self, mydir):
        raise NotImplementedError
