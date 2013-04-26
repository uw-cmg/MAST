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
        if self.is_ready_to_run():
            return
        if self.keywords['program'] == 'vasp':
            self.set_up_vasp_optimize()
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")

        #while not self.is_ready_to_run():
        #    print 'writing files...'
        #    time.sleep(CHECKPERIOD)
            

    def set_up_vasp_incar_dict(self, rep_structure, rep_potcar):
        myd=dict()
        for key, value in self.keywords['program_keys'].iteritems():
            myd[key.upper()]=value
        #myd['IBRION']=self.keywords['program_keys']['ibrion']
        #myd['ISIF']=3
        #myd['NSW']=191
        #myd['NPAR']=4 #hardcoded - needs fixing
        #myd['PREC']="Accurate"
        #myd['ISMEAR']=1
        #myd['SIGMA']=0.2
        #myd['ISPIN']=2
        #myd['LCHARG']="False"
        #myd['LWAVE']="False"
        #myd['NSW']=191
        myd['MAGMOM']=str(len(rep_structure.sites)) + "*5"
        myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(rep_potcar)*1.5
        return myd

    
    def set_up_vasp_optimize(self):
        name = self.keywords['name']
        pospath = os.path.join(name, "POSCAR")
        if self.keywords['structure'] == None:
            opt_poscar = Poscar.from_file(pospath) 
            #parent should have given a structure
        else: #this is an originating run; mast should give it a structure
            opt_poscar = Poscar(self.keywords['structure'])
            self.lock_directory()
            opt_poscar.write_file(pospath)
            self.unlock_directory()
        topkpoints = Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        self.lock_directory()
        topkpoints.write_file(name + "/KPOINTS")
        self.unlock_directory()
        toppotcar = Potcar(symbols=opt_poscar.site_symbols, functional='PBE', sym_potcar_map=None)
        self.lock_directory()
        toppotcar.write_file(name + "/POTCAR")
        self.unlock_directory()
        incar_dict = self.set_up_vasp_incar_dict(opt_poscar.structure, toppotcar)
        topincar = Incar(incar_dict)
        self.lock_directory()
        topincar.write_file(name + "/INCAR")
        self.unlock_directory()
        self.write_submit_script() #MASTFile types already lock/unlock
        return
    
    def write_submit_script(self):
        return BaseIngredient.write_submit_script(self)

    def is_ready_to_run(self):
        return BaseIngredient.is_ready_to_run(self) 

    def run(self, mode='noqsub', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)
    
    # hw 04/15/13 This will be used by scheduler
    def getpath(self):
        return BaseIngredient.getpath(self)
        
