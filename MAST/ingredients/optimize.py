import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import MASTObj
from MAST.ingredients.baseingredient import BaseIngredient
from MAST.utility import MASTError
from MAST.utility import dirutil

import os
import shutil
import subprocess
import time
#TA


class Optimize(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'child_dict': (dict, dict(), 'Dictionary of children'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
    
    def is_complete(self):
        return BaseIngredient.is_complete(self)

    def update_children(self):
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname)
            
    def write_files(self):
        if self.keywords['program'] == 'vasp':
            self.set_up_program_input()
            self.write_submit_script()
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")

        #while not self.is_ready_to_run():
        #    print 'writing files...'
        #    time.sleep(CHECKPERIOD)
            
    def run(self, mode='serial', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)
