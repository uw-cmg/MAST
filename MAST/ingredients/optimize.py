import numpy as np

from pymatgen.core.structure import Structure
#from pymatgen.core.structure_modifier import StructureEditor
#from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar

from MAST.utility import MASTObj
from MAST.ingredients import BaseIngredient
from MAST.utility import MASTError

import os
import shutil
#TA


class Optimize(BaseIngredient):
    def __init__(self, **kwargs):
        allowed_keys = {
                'name' : (str, str(), 'Name of optimization directory'),
                'structure': (Structure, None, 'Pymatgen Structure object'),
                'program': (str, str(), 'DFT program, e.g. "vasp"'),
                'child_dict': (dict, dict(), 'Dictionary of children')
                }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

    def is_complete(self):
        return BaseIngredient.is_complete() #instead call base ingredient complete check

    def update_children(self):
        pass 

    def write_files(self):
        if self.keywords['program'] == 'vasp':
            self.set_up_vasp_optimize()
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")

    def set_up_vasp_incar_dict(self, rep_structure, rep_potcar):
        myd=dict()
        myd['IBRION']=self.keywords['program_keys']['ibrion']
        myd['ISIF']=3
        myd['NSW']=191
        myd['NPAR']=4 #hardcoded - needs fixing
        myd['PREC']="Accurate"
        myd['ISMEAR']=1
        myd['SIGMA']=0.2
        myd['ISPIN']=2
        myd['LCHARG']="False"
        myd['LWAVE']="False"
        myd['NSW']=191
        myd['MAGMOM']=5*len(rep_structure.sites)
        myd['ENCUT']=MAST.ingredients.pmgextend.vasp_extend.get_max_enmax_from_potcar(rep_potcar)*1.5
        return myd

    
    def set_up_vasp_optimize(self):
        name = self.keywords['name']
        if not (self.keywords['structure'] == None):
            pospath = os.path.join(name, "POSCAR")
            opt_poscar = Poscar(self.keywords['structure'])
            opt_poscar.write_file(pospath)
        else:
            opt_poscar = Poscar.from_file(name, "POSCAR")
        topkpoints = pymatgen.io.vaspio.Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        topkpoints.write_file(name + "/KPOINTS")
        toppotcar = pymatgen.io.vaspio.Potcar(symbols=opt_poscar.site_symbols, functional='PAW-GGA', sym_potcar_map=None)
        toppotcar.write_file(name + "/POTCAR")
        topincar = pymatgen.io.vaspio.Incar()
        incar_dict = self.set_up_vasp_incar_dict(opt_poscar.structure, toppotcar)
        topincar.from_dict(incar_dict)
        topincar.write_file(name + "/INCAR")
        return

