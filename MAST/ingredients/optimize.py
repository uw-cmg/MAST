import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar

from MAST.utility import MASTObj
from MAST.ingredients import BaseIngredient

import os
import shutil
import pdb



class Optimize(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
                'dir_name' : (str, str(), 'Name of optimization directory'),
                'ibrion': (str, str(), 'IBRION for optimization (default 2)'),
                'structure_init': (Structure, None, 'Pymatgen Structure object'),
                #'structure_final': (Structure, None, 'Pymatgen Structure object'),#TA never need for opt, right?
                'parent_init': (str, str(), 'Directory of initial state'),
                #'parent_final': (str, str(), 'Directory of final state'), #TA also don't need
                'program': (str, str(), 'DFT program, e.g. "vasp"')
                }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

    def is_complete(self):
        if self.keywords['program'] == "vasp":
            try:
                myoutcar = Outcar(os.path.join(dir_name, "OUTCAR"))
            except (IOError):
                return False
                

            if not myoutcar.run_stats['User time (sec)'] > 0:
                return False
            else:
                return True
        else:
            print "Program " + self.keywords['program'] + " not yet supported for optimize."
            return None
    
    def get_structure_from_parent(self, parentpath):
        if self.keywords['program'] == "vasp":
            structpath = os.path.join(parentpath, "CONTCAR")
            parent_structure = Poscar.from_file(structpath, False).structure
            return parent_structure
        else:
            print "Could not get structure."
            return None

    def get_energy_file_from_parent(self, parentpath, program=""): #TA I think program param should default to using keywrd
        if program == "":
            program = self.keywords['program']
        if program == "vasp":
            return os.path.join(parentpath, "OSZICAR")
        else:
            print "Could not get energy file."
            return None

    def generate_files(self):
        if not self.keywords['parent_init'] == None: #TA could be empty for optimize
            struct_init = self.get_structure_from_parent(self.keywords['parent_init'])
        else: #TA if empty then use the keyword structure
            struct_init = self.keywords['structure_init']
        self.set_up_optimize(struct_init)
        return
    
    def set_up_optimize(self, struct_init):
        dir_name = self.keywords['dir_name']
        if os.path.exists(dir_name):
            print "Directory exists."
            return
        if self.keywords['program'] == 'vasp':
            os.makedir(dir_name)
            opt_poscar = Poscar(struct_init)
            pospath = os.path.join(dir_name, "POSCAR")
            opt_poscar.write_file(pospath)
        else:
            print "Program " + self.keywords['program'] + " not supported for optimize."
            return None

