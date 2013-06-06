
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
from MAST.utility.mastfile import MASTFile
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
        return BaseIngredient.is_complete(self)

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
            BaseIngredient.set_up_program_input_neb(self, image_structures)
            self.place_parent_energy_files()
            self.write_submit_script()
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
        numstr = numstr[:numstr.find("_")]
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
        structure_list = struct_init.interpolate(struct_fin, int(self.keywords['program_keys']['images']) + 1 )
        return structure_list

    def place_parent_energy_files(self):
        """Assume parents have written files parent_energy_<N>.
            Copy these files into the 00 and 0N directories.
        """
        numparents = self.get_my_numbers()
        if numparents == None:
            raise MASTError(self.__class__.__name__,"No parents given!")
        pfpath1=os.path.join(self.keywords['name'], "parent_energy_" + numparents[0])
        pfpath2=os.path.join(self.keywords['name'], "parent_energy_" + numparents[1])
        pffile1=MASTFile(pfpath1)
        pffile2=MASTFile(pfpath2)
        pffile1.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self,1)) #MASTFile contains directory locking.
        pffile2.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, 2))
        return

    def run(self, mode='serial', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)

