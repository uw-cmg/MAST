
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
        for childname in sorted(self.keywords['child_dict'].iterkeys()):
            impath = os.path.join(self.keywords['name'], str(myct).zfill(2))
            self.forward_parent_structure(impath, childname)       
            myct = myct + 1

    def write_files(self):
        parentstructures = self.get_parent_structures()
        image_structures = self.do_interpolation(parentstructures)
        if image_structures == None:
            raise MASTError(self.__class__.__name__,"Bad number of images")
        if self.keywords['program'] == 'vasp':
            self.set_up_vasp_neb(image_structures)
        else:
            raise MASTError(self.__class__.__name__,"Program not supported. No setup accomplished.")
        return
   
    def get_my_numbers(self):
        """For neb in the format <sys>_neb<N>-<N> return the two
            key numbers identifying the defect sites.
        """
        tempname = self.keywords['name'].lower()
        numstr = tempname[tempname.find("neb")+3:]
        if numstr.find("neb") == -1:
            pass
        else:
            numstr = numstr[numstr.find("neb")+3:] #allow 'neb' to be in <sys>
        numparents = numstr.split('-')
        if not(len(numparents) == 2):
            raise MASTError(self.__class__.__name__,"Error: number of parents is not two for NEB in " + self.keywords['name'])
            return None
        return numparents


    def get_parent_structures(self):
        """Assume that parents have written two files,
            named 'parent_structure_<N>'. 
            For VASP these are CONTCAR-type files.
            Returns:
                [struct_init, struct_fin]: pymatgen Structure objects
        """
        numparents = self.get_my_numbers()
        pfpath1=os.path.join(self.keywords['name'], "parent_structure_" + numparents[0])
        if not os.path.isfile(pfpath1):
            raise MASTError(self.__class__.__name__,"Error: no parent file at" + pfpath1)
        struct_init = BaseIngredient.get_structure_from_file(self, pfpath1)
        pfpath2=os.path.join(self.keywords['name'], "parent_structure_" + numparents[1])
        if not os.path.isfile(pfpath2):
            raise MASTError(self.__class__.__name__,"Error: no parent file at" + pfpath2)
        struct_fin = BaseIngredient.get_structure_from_file(self, pfpath2)
        return [struct_init, struct_fin]

    def do_interpolation(self, parentstructures):
        """Do interpolation."""
        if parentstructures == None:
            raise MASTError(self.__class__.__name__,"Bad number of parent paths.")
        struct_init = parentstructures[0]
        struct_fin = parentstructures[1]
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
    

    def place_parent_energy_files_vasp(self):
        """Assume parents have written files parent_energy_<N>.
            Copy these files into the 00 and 0N directories.
        """
        numparents = self.get_my_numbers()
        if numparents == None:
            return
        pfpath1=os.path.join(self.keywords['name'], "parent_energy_" + numparents[0])
        if not os.path.isfile(pfpath1):
            raise MASTError(self.__class__.__name__,"Error: no parent file at" + pfpath1)
        shutil.copy(pfpath1, os.path.join(self.keywords['name'],"00","OSZICAR"))
        pfpath2=os.path.join(self.keywords['name'], "parent_energy_" + numparents[1])
        if not os.path.isfile(pfpath2):
            raise MASTError(self.__class__.__name__,"Error: no parent file at" + pfpath2)
        shutil.copy(pfpath2, os.path.join(self.keywords['name'],str(self.keywords['program_keys']['images']+1).zfill(2),"OSZICAR"))
        return

    def set_up_vasp_folders(self, image_structures):
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
            imct = imct + 1
        self.place_parent_energy_files_vasp()
        return
        

    def set_up_vasp_neb(self, image_structures):
        self.set_up_vasp_folders(image_structures)
        dir_name = self.keywords['name']
        topkpoints = Kpoints.monkhorst_automatic(kpts=(4,4,4),shift=(0,0,0))
        topkpoints.write_file(dir_name + "/KPOINTS")
        toppotcar = Potcar(symbols=Poscar(image_structures[0]).site_symbols, functional='PBE', sym_potcar_map=None)
        toppotcar.write_file(dir_name + "/POTCAR")
        incar_dict = self.set_up_vasp_incar_dict(image_structures[0], toppotcar)
        topincar = Incar(incar_dict)
        topincar.write_file(dir_name + "/INCAR")
        return


