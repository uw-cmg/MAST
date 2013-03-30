import numpy as np

import pymatgen
from pymatgen.core.structure import Structure
from pymatgen.io.vaspio import Poscar
from pymatgen.io.vaspio import Outcar

from MAST.ingredients.pmgextend import vasp_extensions
from MAST.utility import MASTObj
from MAST.ingredients import BaseIngredient

import os
import shutil
#import pdb
#TTM

class PerformNEB(BaseIngredient):
    def __init__(self, **kwargs):
        #pdb.set_trace()
        allowed_keys = {
                'dir_name' : (str, str(), 'Name of NEB directory'),
                'parent_init': (str, str(), 'Directory of initial state'),
                'parent_final': (str, str(), 'Directory of final state'),
                'images': (int, 0, 'Number of images'),
                'program': (str, str(), 'DFT program, e.g. "vasp"')
                }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)
        #prog_kwarg_dict = options.get_item('vasp',ingredient_name)

    def is_complete(self):
        return BaseIngredient.images_complete(self)

    def generate_files(self):
        image_structures = self.do_interpolation()
        if self.keywords['program'] == 'vasp':
            self.set_up_vasp_neb(image_structures)
        else:
            print "Program not supported. No setup accomplished."
            return
        return
   
    def do_interpolation(self):
        struct_init = None
        struct_fin = None
        struct_init = BaseIngredient.get_structure_from_parent(self, self.keywords['parent_init'])
        struct_fin = BaseIngredient.get_structure_from_parent(self, self.keywords['parent_final'])
        if (struct_init == None) or (struct_fin == None):
            print "Error getting initial or final parent structure."
            return
        structure_list = struct_init.interpolate(struct_fin, self.keywords['images'] + 1 )
        return structure_list

    def set_up_vasp_incar_dict(self, rep_structure, rep_potcar):

        myd=dict()
        myd['IBRION']=1
        myd['POTIM']=0.5
        myd['ISIF']=2
        myd['LCLIMB']=True
        myd['NSW']=191
        myd['NPAR']=4
        myd['PREC']="Accurate"
        myd['ISMEAR']=0
        myd['SIGMA']=0.05
        myd['ISPIN']=2
        myd['IMAGES']=self.keywords['images']
        myd['LCHARG']="False"
        myd['LWAVE']="False"
        myd['NSW']=191
        myd['MAGMOM']=str(len(rep_structure.sites)) + "*5"
        myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(rep_potcar)*1.5
        return myd

    def set_up_vasp_folders(self, image_structures):
        dir_name = self.keywords['dir_name']
        if os.path.exists(dir_name):
            print "Directory exists."
            return
        os.makedirs(dir_name)
        imct=0
        while imct <= self.keywords['images'] + 1:
            imposcar = Poscar(image_structures[imct])
            num_str = str(imct).zfill(2)
            impath = os.path.join(dir_name, num_str)
            impospath = os.path.join(dir_name, "POSCAR_" + num_str)
            imposcar.write_file(impospath)
            os.makedirs(impath)
            imposcar.write_file(os.path.join(impath, "POSCAR"))
            if imct == 0:
                shutil.copy(os.path.join(self.keywords['parent_init'],"OSZICAR"),impath)
            elif imct == self.keywords['images'] + 1:
                shutil.copy(os.path.join(self.keywords['parent_final'],"OSZICAR"),impath)
            imct = imct + 1
        return
        

    def set_up_vasp_neb(self, image_structures):
        self.set_up_vasp_folders(image_structures)
        dir_name = self.keywords['dir_name']
        topkpoints = pymatgen.io.vaspio.Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        topkpoints.write_file(dir_name + "/KPOINTS")
        toppotcar = pymatgen.io.vaspio.Potcar(symbols=Poscar(image_structures[0]).site_symbols, functional='PBE', sym_potcar_map=None)
        toppotcar.write_file(dir_name + "/POTCAR")
        incar_dict = self.set_up_vasp_incar_dict(image_structures[0], toppotcar)
        topincar = pymatgen.io.vaspio.Incar(incar_dict)
        topincar.write_file(dir_name + "/INCAR")
        return


