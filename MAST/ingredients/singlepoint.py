import numpy as np

from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar
from pymatgen.io.vaspio import Kpoints
from pymatgen.io.vaspio import Potcar
from pymatgen.io.vaspio import Incar

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import MASTObj
from MAST.ingredients import BaseIngredient
from MAST.utility import MASTError
from MAST.utility import dirutil

import os
import shutil
import subprocess
import time
#TA


class SinglePoint(BaseIngredient):
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
# Redundant, keep it for now
        return BaseIngredient.is_complete(self) #instead call base ingredient complete check

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

        while not self.is_ready_to_run():
            print 'writing files...'
            time.sleep(CHECKPERIOD)

    def set_up_vasp_incar_dict(self, rep_structure, rep_potcar):
        myd=dict()
        for key, value in self.keywords['program_keys'].iteritems():
            myd[key.upper()]=value
        myd['MAGMOM']=str(len(rep_structure.sites)) + "*5"
        myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(rep_potcar)*1.5
        return myd

    def set_up_vasp_singlepoint(self):
        name = self.keywords['name']
        pospath = os.path.join(name, "POSCAR")
        if self.keywords['structure'] == None:
            opt_poscar = Poscar.from_file(pospath) 
            #parent should have given a structure
        else: #this is an originating run; mast should give it a structure
            opt_poscar = Poscar(self.keywords['structure'])
            opt_poscar.write_file(pospath)
        topkpoints = Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        topkpoints.write_file(name + "/KPOINTS")
        toppotcar = Potcar(symbols=opt_poscar.site_symbols, functional='PBE', sym_potcar_map=None)
        toppotcar.write_file(name + "/POTCAR")
        incar_dict = self.set_up_vasp_incar_dict(opt_poscar.structure, toppotcar)
        topincar = Incar(incar_dict)
        topincar.write_file(name + "/INCAR")
        self.write_submit_script()
        return

    def write_submit_script(self):
        myfile = open(dirutil.get_mast_install_path() + '/submit/node1.sh','rb')
        mylines=myfile.readlines()
        myfile.close()
        bname = os.path.basename(self.keywords['name'])
        myct=0
        while myct < len(mylines):
            if "#PBS -N" in mylines[myct]:
                mylines[myct] = "#PBS -N " + bname + '\n'
            myct=myct+1
        mywrite = open(self.keywords['name']+'/submit.sh','wb')
        mywrite.writelines(mylines)
        mywrite.close()
        return

    def is_ready_to_run(self):
        name = self.keywords['name']
        notready=0
        if not(os.path.isfile(name + "/KPOINTS")):
            notready = notready + 1
        if not(os.path.isfile(name + "/POTCAR")):
            notready = notready + 1
        if not(os.path.isfile(name + "/INCAR")):
            notready = notready + 1
        if not(os.path.isfile(name + "/POSCAR")):
            notready = notready + 1
        if not(os.path.isfile(name + "/submit.sh")):
            notready = notready + 1
        if notready > 0:
            return False
        else:
            return True

    def run(self, mode='noqsub', curdir=os.getcwd()):
        if not (self.is_ready_to_run()):
            # we need a MAST Warning class
            raise
 
        if mode is 'noqsub':
            curdir = os.getcwd()
            os.chdir(self.keywords['name'])
            programpath = '//share/apps/vasp5.2_cNEB' # This should be replaced by more general way.
            p = subprocess.call([programpath])
            os.chdir(curdir)
            
        elif mode is 'serial':
            curdir = os.getcwd()
            os.chdir(self.keywords['name'])
            runme = subprocess.Popen('qsub submit.sh', shell=True, stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)
            os.chdir(curdir)
            # for scheduling other jobs
            #runme.wait()
            return
    
    # hw 04/15/13 This will be used by scheduler
    def getpath(self):
        '''getpath returns the directory of the ingredient'''
        return self.keywords['name']

