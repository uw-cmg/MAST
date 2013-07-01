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

class NEB(BaseIngredient):
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
        parentimagestructures = self.get_parent_image_structures()
        if len(parentimagestructures) == 0:
            image_structures = self.do_interpolation(parentstructures)
        else:
            image_structures = list()
            image_structures.append(parentstructures[0])
            image_structures.extend(parentimagestructures)
            image_structures.append(parentstructures[1])
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
        """For neb in the format <sys>_neb_<N>-<N>_... return the two
            labels identifying the defect sites.
        """
        myname = os.path.basename(self.keywords['name'])
        tempname = myname.lower()
        lsplit = tempname.split('_') #last underscore segment
        nebidx = lsplit.index('neb')
        lseg = lsplit[nebidx + 1] #Follows after "neb"
        llist = lseg.split('-')
        if not(len(llist) == 2):
            raise MASTError(self.__class__.__name__, 
                "Error: number of NEB parents <> 2 in " + myname+":"+str(llist))
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
        sortedstruc_base = self.keywords['structure'].get_sorted_structure()
        struct_ed = StructureEditor(sortedstruc)
        struct_ed_base = StructureEditor(sortedstruc_base)
        # We are not actually translating the sites; just get them into unit
        # cell
        #TTM DEBUG REMOVE THIS struct_ed.translate_sites(range(0,len(sortedstruc.sites)),np.zeros(3),True)
        nebidx = list()

        elemstarts = self.get_element_indices(sortedstruc)
        mylabels = self.get_my_labels()
        mylabel = '-'.join(mylabels)
        #print "TTM DEBUG: my neb lines", self.neblines[mylabel]
        for nebline in self.neblines[mylabel]:
            usebase=0
            nebdict = self.parse_neb_line(nebline)
            mycoord = nebdict['coord'][initorfin]
            index = find_in_coord_list(sortedstruc.frac_coords, mycoord, atol)
            print 'TTM DEBUG: index: ', index
            if len(index) == 0: #try the base structure, for vacancies
                usebase=1
                index = find_in_coord_list(sortedstruc_base.frac_coords,
                    mycoord, atol)
            if len(index) == 0:
                raise MASTError(self.__class__.__name__, "No coordinate found matching %s" % mycoord)
            if usebase==0:
                nebidx.append(index[0]) #only take first site?
                mysite = sortedstruc.sites[index[0]]
                myelem = MAST.data.atomic_number[mysite.species_string]
                struct_ed.delete_site(index)
                struct_ed.insert_site(elemstarts[myelem], mysite.specie,
                                        mysite.frac_coords)
            else:
                nebidx.append(index[0]) #only take first site?
                mysite = sortedstruc_base.sites[index[0]]
                myelem = MAST.data.atomic_number[mysite.species_string]
                struct_ed_base.delete_site(index)
                struct_ed_base.insert_site(elemstarts[myelem], mysite.specie,
                                        mysite.frac_coords)

        if not len(nebidx) == len(self.neblines[mylabel]):
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
        """Parse an NEB atomic movement line which has the following format:
                    Cr, 0 0 0, 0.5 0.5 0.5 in the input file, and
                    ['Cr','0 0 0','0.5 0.5 0.5'] here as nebline
                    element, init_coord, fin_coord
            Args:
                line <list of str>: NEB atomic movement line
            Returns:
                linedict <dict>: 
                                 ['element'] <int> = element's atomic number
                                 ['coord'][0] <np.array> = initial coordinates
                                 ['coord'][1] <np.array> = final coordinates
        """
        nebdict=dict()
        #print "TTM DEBUG: nebline", nebline
        nebelem = str(nebline[0].strip()).title()
        import MAST.data
        nebdict['element'] = MAST.data.atomic_number[nebelem]
        nebdict['coord'] = dict()
        nebdict['coord'][0] = np.array(nebline[1].split(), dtype='float')
        nebdict['coord'][1] = np.array(nebline[2].split(), dtype='float')
        return nebdict
        
    def get_parent_structures(self):
        """Assume that parents have written files
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

    def get_parent_image_structures(self):
        """A low-mesh NEB may have written files
            named 'parent_structure_<N-N>_0N'. 
            For VASP these are CONTCAR-type files.
            Returns:
                list of <Structure>: list of pymatgen Structure objects
        """
        header = "parent_structure_"
        numim = self.keywords['program_keys']['images']
        imct = 1
        imstrs=list()
        while imct <= numim:
            pfpath=""
            for myfile in os.listdir(self.keywords['name']):
                if (header in myfile) and (str(imct).zfill(2) in myfile):
                    pfpath = os.path.join(self.keywords['name'],myfile)
            if pfpath == "":
                pass
            else:
                struct_im = BaseIngredient.get_structure_from_file(self, pfpath)
                sorted_im = self.sort_structure_and_neb_lines(struct_im, 0) 
                imstrs.append(sorted_im)
            imct = imct + 1
        if len(imstrs) > 0 and not (len(imstrs) == numim):
            raise MASTError(self.__class__.__name__, "Incomplete number of forwared images found!")
        return imstrs
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

