import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar

from MAST.utility.mastobj import MASTObj
from MAST.ingredients.libingredients import BaseIngredient

import os
import shutil
import pdb



class PerformNEB(MASTObj):
    def __init__(self, **kwargs):
        #pdb.set_trace()
        allowed_keys = {
                'dir_name' : (str, str(), 'Name of NEB directory'),
                'ibrion': (str, str(), 'IBRION for NEB (default 1)'),
                'potim': (str, str(), 'POTIM for NEB (default 0.5)'),
                'climbing': (str, str(), 'Climbing NEB? T/F'),
                'images': (int, 4, 'Number of interpolated images + 1'),
                'structure_init': (Structure, None, 'Pymatgen Structure object'),
                'structure_final': (Structure, None, 'Pymatgen Structure object'),
                'parent_init': (str, str(), 'Directory of initial state'),
                'parent_final': (str, str(), 'Directory of final state'),
                'program': (str, str(), 'DFT program, e.g. "vasp"')
                }
        MASTObj.__init__(self, allowed_keys, **kwargs)
        #BaseIngredient.__init__(self,**kwargs)

    def is_complete(self):
        if self.keywords['program'] == "vasp":
            imct=1
            while imct < self.keywords['images']:
                num_str = str(imct).zfill(2)
                impath = os.path.join(self.keywords['dir_name'], num_str)
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
        else:
            print "Program " + self.keywords['program'] + " not supported for NEB."
            return None
    
    def get_structure_from_parent(self, parentpath):
        if self.keywords['program'] == "vasp":
            structpath = os.path.join(parentpath, "CONTCAR")
            parent_structure = Poscar.from_file(structpath, False).structure
            return parent_structure
        else:
            print "Could not get structure."
            return None

    def get_energy_file_from_parent(self, parentpath, program):
        if program == "vasp":
            return os.path.join(parentpath, "OSZICAR")
        else:
            print "Could not get energy file."
            return None

    def generate_files(self):
        struct_init = self.get_structure_from_parent(self.keywords['parent_init'])
        struct_fin = self.get_structure_from_parent(self.keywords['parent_final'])
        image_structures = struct_init.interpolate(struct_fin, self.keywords['images'])
        self.set_up_neb(image_structures)
        return
    
    def set_up_neb(self, image_structures):
        dir_name = self.keywords['dir_name']
        if os.path.exists(dir_name):
            print "Directory exists."
            return
        if self.keywords['program'] == 'vasp':
            os.makedirs(dir_name)
            imct=0
            while imct <= self.keywords['images']:
                imposcar = Poscar(image_structures[imct])
                num_str = str(imct).zfill(2)
                impath = os.path.join(dir_name, num_str)
                impospath = os.path.join(dir_name, "POSCAR_" + num_str)
                imposcar.write_file(impospath)
                os.makedirs(impath)
                imposcar.write_file(os.path.join(impath, "POSCAR"))
                if imct == 0:
                    shutil.copy(os.path.join(self.keywords['parent_init'],"OSZICAR"),impath)
                elif imct == self.keywords['images']:
                    shutil.copy(os.path.join(self.keywords['parent_final'],"OSZICAR"),impath)
                imct = imct + 1
        else:
            print "Program " + self.keywords['program'] + " not supported for NEB."
            return None


