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
from MAST.utility.mastfile import MASTFile
import os
import shutil
#import pdb
#TTM

class PerformNEB(BaseIngredient):
    """
        Attributes:
        self.labels <list of str>: list of labels labeling the NEB.
                                    self.labels[0] is the initial state label
                                    self.labels[1] is the final state label
        self.neblines <list of str>: list of NEB atom motion lines for this NEB
    """
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
        self.neblines = self.keywords['program_keys']['neblines']
        self.labels = self.get_my_labels()

    def is_complete(self):
        return BaseIngredient.is_complete(self)

    def update_children(self):
        """Update children by forwarding parent structure into image folders.
        """
        myct=1
        impath=""
        for childname in sorted(self.keywords['child_dict'].iterkeys()):
            impath = os.path.join(self.keywords['name'], str(myct).zfill(2))
            self.forward_parent_structure(impath, childname)       
            myct = myct + 1

    def write_files(self):
        """Get the parent structures, sort and match atoms, and interpolate.
            Write images to the appropriate folders.
        """
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
   
    def get_my_labels(self):
        """For neb in the format <sys>_neb_<N>-<N> return the two
            labels identifying the defect sites.
        """
        myname = self.keywords['name']
        tempname = myname.lower()
        lsplit = tempname.split('_')[-1] #last underscore segment
        llist = lsplit.split('-')
        if not(len(llist) == 2):
            raise MASTError(self.__class__.__name__, 
                "Error: number of NEB parents <> 2 in " + myname)
            return None
        return llist


    def sort_structure_and_neb_lines(self, mystruc, initorfin):
        """Sort the structure in an approved way:
            1. Edge-coorect the fractional structure so for example 
                -0.05 and 0.95 are the same.
            2. Remove lines which closely match the neb moving lines in the 
            NEB section.
            3. Sort the structure by element and coordinate.
            4. Prepend the NEB moving lines to their element sections.
            5. Return the modified structure.

            mystruc <Structure>   : pymatgen structure
            initorfin <int>       : 0 = initial cfg, 1 = final cfg
        """
        from pymatgen.core.structure_modifier import StructureEditor
        from pymatgen.util.coord_utils import find_in_coord_list
        import MAST.data
        atol = 0.1 # Need fairly large tolerance to account for relaxation.
        sortedstruc = mystruc.get_sorted_structure()
        struct_ed = StructureEditor(sortedstruc)
        # We are not actually translating the sites; just get them into unit
        # cell
        struct_ed.translate_sites(range(0,len(sortedstruc.sites)),np.zeros(3),True)
        nebidx = list()

        elemstarts = self.get_element_indices(sortedstruc)

        for nebline in self.neblines:
            nebdict = self.parse_neb_line(nebline)
            index = find_in_coord_list(sortedstruc.frac_coords,
                    nebdict['coord'][initorfin],atol)
            nebidx.append(index)
            mysite = sortedstruc.sites[index]
            myelem = MAST.data.atomic_number[mysite.species_string]
            struct_ed.delete_site(index)
            struct_ed.insert_site(elemstarts[myelem], mysite.specie,
                                    mysite.frac_coords)
        if not len(nebidx) == len(self.neblines):
            raise MASTError(self.__class__.__name__, "Not all NEB lines found.")
        return struct_ed.modified_structure        

    def get_element_indices(self, sortedstruc):
        """From a sorted structure, get the element indices
            Args:
                sortedstruc <Structure>: pymatgen structure sorted by
                                            electronegativity
            Returns:
                elstart <dict>: element index dictionary, in the format:
                            elstart[element number] = starting index, e.g.
                            elstart[24] = 0
        """
        elstart = dict()
        numlist = sortedstruc.atomic_numbers
        listlen = len(numlist)
        idx = 0
        while idx < listlen:
            atomnum = numlist[idx]
            if not atomnum in elstart.keys():
                elstart[atomnum] = idx
            idx = idx + 1
        return elstart

    def parse_neb_line(self, nebline):
        """Parse an NEB line with the following format:
                    1-2, Cr, 0 0 0, 0.5 0.5 0.5
                    label, element, init_coord, fin_coord
            Args:
                line <str>: NEB atomic movement line
            Returns:
                linedict <dict>: ['label'] <str> = label
                                 ['element'] <int> = element's atomic number
                                 ['coord'][0] <np.array> = initial coordinates
                                 ['coord'][1] <np.array> = final coordinates
        """
        nebsplit = nebline.split(',')
        nebdict=dict()
        nebdict['label'] = str(nebsplit[0].strip())
        nebelem = str(nebsplit[1].strip())
        import MAST.data
        nebdict['element'] = MAST.data.atomic_number[nebelem]
        nebdict['coord'] = dict()
        nebdict['coord'][0] = np.array(nebsplit[2].split(), dtype='float')
        nebdict['coord'][1] = np.array(nebsplit[3].split(), dtype='float')
        return nebdict
        
    def get_parent_structures(self):
        """Assume that parents have written two files,
            named 'parent_structure_<N>'. 
            For VASP these are CONTCAR-type files.
            Returns:
                [struct_init, struct_fin]: pymatgen Structure objects
        """
        header = os.path.join(self.keywords['name'], "parent_structure_")
        pfpath_init = header + self.labels[0]
        pfpath_fin = header + self.labels[1]
        if not os.path.isfile(pfpath_init):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_init)
        if not os.path.isfile(pfpath_fin):
            raise MASTError(self.__class__.__name__,
                "Error: no parent file at" + pfpath_fin)
        struct_init = BaseIngredient.get_structure_from_file(self, pfpath_init)
        struct_fin = BaseIngredient.get_structure_from_file(self, pfpath_fin)
        sorted_init = self.sort_structure_and_neb_lines(struct_init, 0) 
        sorted_fin = self.sort_structure_and_neb_lines(struct_fin, 1)
        return [sorted_init, sorted_fin]

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
        header = os.path.join(self.keywords['name'], "parent_energy_")
        pfpath1= header + self.labels[0]
        pfpath2= header + self.labels[1]
        pffile1=MASTFile(pfpath1)
        pffile2=MASTFile(pfpath2)
        pffile1.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self,1)) #MASTFile contains directory locking.
        pffile2.to_file(BaseIngredient.get_path_to_write_neb_parent_energy(self, 2))
        return

    def run(self, mode='serial', curdir=os.getcwd()):
        if not (self.is_ready_to_run()): #This check must occur here in case is_ready_to_run is ever overridden directly in the class.
            raise MASTError(self.__class__.__name__, "Asked to run job before job was ready.")
        return BaseIngredient.run(self, mode)

