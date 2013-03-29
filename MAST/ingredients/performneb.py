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
                'name' : (str, str(), 'Name of NEB directory'),
                'program': (str, str(), 'DFT program, e.g. "vasp"'),
                'program_keys': (dict, dict(), 'Program-specific keywords'),
                'child_dict':(dict, dict(), 'Dictionary of child names'),
                'structure': (Structure, None, 'Pymatgen-type structure')
                }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

    def is_complete(self):
        return BaseIngredient.images_complete(self)

    def update_children(self):
        myct=1
        impath=""
        for childname in child_dict:
            impath = os.path.join(self.keywords['name'], str(myct).zfill(2))
            self.forward_parent_structure(impath, childname)       
            myct = myct + 1

    def write_files(self):
        parentpaths = self.get_parent_paths()
        image_structures = self.do_interpolation(parentpaths)
        if image_structures == None:
            print "Bad number of images."
            return None
        if self.keywords['program'] == 'vasp':
            self.set_up_vasp_neb(image_structures, parentpaths)
        else:
            print "Program not supported. No setup accomplished."
            return
        return
   

    def get_parent_paths(self):
        """Assumes the NEB is named <sys>_NEB<N>-<N> and that the
            parents have written their pathnames into parent_path_<N>.
        """
        tempname = self.keywords['name'].lower()
        numstr = tempname[tempname.find("neb")+3:]
        if numstr.find("neb") == -1:
            pass
        else:
            numstr = numstr[numstr.find("neb")+3:] #allow 'neb' to be in sys name
        numpaths = numstr.split('-')
        myfiles=os.listdir(self.keywords['name'])
        parentpaths=[]
        pfile=""
        print numpaths
        for onenum in numpaths:
            pfpath=os.path.join(self.keywords['name'], "parent_path_" + onenum)
            if os.path.isfile(pfpath):
                pfile = open(pfpath, 'rb')
                parentpaths.append(pfile.readlines()[0].strip())
                pfile.close()
        print parentpaths
        return parentpaths

    def do_interpolation(self, parentpaths):
        struct_init = None
        struct_fin = None
        if not (len(parentpaths) == 2):
            print "Bad number of parent paths."
            return 
        struct_init = BaseIngredient.get_structure_from_parent(self, parentpaths[0])
        struct_fin = BaseIngredient.get_structure_from_parent(self, parentpaths[1])
        if (struct_init == None) or (struct_fin == None):
            print "Error getting initial or final parent structure."
            return
        structure_list = struct_init.interpolate(struct_fin, self.keywords['program_keys']['images'] + 1 )
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
        myd['IMAGES']=self.keywords['program_keys']['images']
        myd['LCHARG']="False"
        myd['LWAVE']="False"
        myd['NSW']=191
        myd['MAGMOM']=str(len(rep_structure.sites)) + "*5"
        myd['ENCUT']=vasp_extensions.get_max_enmax_from_potcar(rep_potcar)*1.5
        return myd

    def set_up_vasp_folders(self, image_structures, parentpaths):
        imct=0
        while imct <= self.keywords['program_keys']['images'] + 1:
            imposcar = Poscar(image_structures[imct])
            num_str = str(imct).zfill(2)
            impath = os.path.join(self.keywords['name'], num_str)
            impospath = os.path.join(self.keywords['name'], "POSCAR_" + num_str)
            imposcar.write_file(impospath)
            try:
                os.makedirs(impath)
            except OSError:
                print "Directory at", impath, "already exists."
                return None
            imposcar.write_file(os.path.join(impath, "POSCAR"))
            if imct == 0:
                shutil.copy(os.path.join(parentpaths[0],"OSZICAR"),impath)
            elif imct == self.keywords['program_keys']['images'] + 1:
                shutil.copy(os.path.join(parentpaths[1],"OSZICAR"),impath)
            imct = imct + 1
        return
        

    def set_up_vasp_neb(self, image_structures, parentpaths):
        self.set_up_vasp_folders(image_structures, parentpaths)
        dir_name = self.keywords['name']
        topkpoints = pymatgen.io.vaspio.Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        topkpoints.write_file(dir_name + "/KPOINTS")
        toppotcar = pymatgen.io.vaspio.Potcar(symbols=Poscar(image_structures[0]).site_symbols, functional='PBE', sym_potcar_map=None)
        toppotcar.write_file(dir_name + "/POTCAR")
        incar_dict = self.set_up_vasp_incar_dict(image_structures[0], toppotcar)
        topincar = pymatgen.io.vaspio.Incar(incar_dict)
        topincar.write_file(dir_name + "/INCAR")
        return


