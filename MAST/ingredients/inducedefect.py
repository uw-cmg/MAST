import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list
from pymatgen.io.vaspio import Poscar
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.ingredients.baseingredient import BaseIngredient
import os

class InduceDefect(BaseIngredient):
    def __init__(self, **kwargs):
        #TTM move allowed_keys into __init__ and match order of other sections
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'child_dict': (dict, dict(), 'Dictionary of children'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }
        BaseIngredient.__init__(self, allowed_keys, **kwargs)

        #if (self.keywords['coordtype'] == 'cartesian'):
        #    self.keywords['position'] = self._cart2frac(self.keywords['position'])

    def induce_defect(self):
        """Creates a defect, and returns the modified structure
            mast_defect is a dictionary like this: 
            'defect1': {'symbol': 'cr', 'type': 'interstitial', 
                        'coordinates': array([ 0. ,  0.5,  0. ])}}
            'defect2': {'symbol': 'o', 'type': 'vacancy', 
                        'coordinates': array([ 0.25,  0.25,  0.25])}
            'coord_type': 'fractional' 
        """
        #TTM change to use mast_defects passed in from inputparser
        print "SELF.KEYWORDS: ",self.keywords
        print "SELF.KEYWORDS PROGRAM_KEYS: ",self.keywords['program_keys']
        myidx = "defect" + str(self.get_my_number()-1) #ex: if you are defect1, index to ['mast_defects']['defect0']
        mydict = self.keywords['program_keys']['mast_defects'][myidx]
        struct_ed = StructureEditor(self.keywords['structure']) #should be updated using get_new_structure)

        if (self.keywords['program_keys']['mast_defects']['coord_type'] == 'cartesian'):
            mydict['coordinates'] = self._cart2frac(mydict['coordinates'])

        if (mydict['type'].lower() == 'vacancy'):
            print 'Making a vacancy!'

            index = find_in_coord_list(self.keywords['structure'].frac_coords,
                                       mydict['coordinates'],
                                       atol=1e-04)

            struct_ed.delete_site(index)
        elif (mydict['type'].lower() == 'interstitial'):
            print 'Making an interstitial!'

            struct_ed.append_site(mydict['symbol'],
                                  mydict['coordinates'],
                                  coords_are_cartesian=False,
                                  validate_proximity=False) # Should be set to True!
        elif (mydict['type'].lower() == 'antisite'):
            print 'Making an antisite!'

            index = find_in_coord_list(self.keywords['structure'].frac_coords,
                                       mydict['coordinates'],
                                       atol=1e-04)

            struct_ed.replace_site(index, mydict['symbol'])
        else:
            raise RuntimeError('Defect type %s not supported' % self.keywords['defecttype'])

        return struct_ed.modified_structure

    def _cart2frac(self, position):
        """Converts between cartesian coordinates and fractional coordinates"""
        return self.keywords['structure'].lattice.get_fractional_coords(self.keywords['position'])

    #def write_input(self):
    #    """Writes the defected geometry to a file"""

    #def _make_directory(self, directory):
    #    if os.path.exists(directory):
    #        print 'Directory %s exists!' % directory
    #    else:
    #        os.path.makedirs(directory)
    
    def write_files(self):
        name=self.keywords['name']
        self.get_new_structure()
        modified_structure = self.induce_defect()
        if self.keywords['program'] == 'vasp':
            myposcar = Poscar(modified_structure)
            self.lock_directory()
            myposcar.write_file(name + '/CONTCAR')
            self.unlock_directory()
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")
        return
    
    def is_ready_to_run(self):
        return True

    def run(self, mode='noqsub'):
        return True

    def is_complete(self):
        if self.directory_is_locked():
            return False
        if self.keywords['program'] == 'vasp':
            if os.path.exists(self.keywords['name'] +'/CONTCAR'):
                return True
            else:
                return False
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")


    def update_children(self):
        for childname in self.keywords['child_dict'].iterkeys():
            self.forward_parent_structure(self.keywords['name'], childname)

    def get_new_structure(self):
        if self.keywords['program'] == 'vasp':
            myposcar = Poscar.from_file(self.keywords['name'] + "/POSCAR")
            self.keywords['structure'] = myposcar.structure
        else:
            raise MASTError(self.__class__.__name__, "Program not supported.")
        return
    
    def get_my_number(self):
        """For defect in the format <sys>_induce_defect<N>, return N.
        """
        tempname = self.keywords['name'].lower()
        numstr = tempname[tempname.find("defect")+6:]
        if numstr.find("defect") == -1:
            pass
        else:
            numstr = numstr[numstr.find("defect")+6:] #allow 'defect' to be in <sys>
        defectnum = int(numstr)
        return defectnum

