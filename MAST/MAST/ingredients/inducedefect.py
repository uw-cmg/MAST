import numpy as np

from pymatgen.core.sites import PeriodicSite
from pymatgen.core.structure import Structure
from pymatgen.core.structure_modifier import StructureEditor
from pymatgen.util.coord_utils import find_in_coord_list

from MAST.utility.mastobj import MASTObj

allowed_keys = {'atom': (str, str(), 'Atom for the defect'),
                'position': (tuple, tuple(), 'Position of the defect'),
                'coordtype': (str, 'fractional', 'Coordinate type (cartesian or fractional'),
                'defecttype': (str, str(), 'Intersitial or vacancy'),
                'structure': (Structure, None, 'Pymatgen Structure object'),
               }

class InduceDefect(MASTObj):
    def __init__(self, **kwargs):
        MASTObj.__init__(self, allowed_keys, **kwargs)

        if (self.keywords['coordtype'] == 'cartesian'):
            self.keywords['position'] = self._cart2frac(self.keywords['position'])

    def induce_defect(self):
        """Creates a defect, and returns the modified structure
        """
        struct_ed = StructureEditor(self.keywords['structure'])

        if (self.keywords['coordtype'] == 'cartesian'):
            self.keywords['position'] = self._cart2frac(self.keywords['position'])

        if (self.keywords['defecttype'].lower() == 'vacancy'):
            print 'Making a vacancy!'

            index = find_in_coord_list(self.keywords['structure'].frac_coords,
                                       self.keywords['position'],
                                       atol=1e-04)

            struct_ed.delete_site(index)
        elif (self.keywords['defecttype'].lower() == 'interstitial'):
            print 'Making an interstitial!'

            struct_ed.append_site(self.keywords['atom'],
                                  self.keywords['position'],
                                  coords_are_cartesian=False,
                                  validate_proximity=False) # Should be set to True!
        elif (self.keywords['defecttype'].lower() == 'antisite'):
            print 'Making an antisite!'

            index = find_in_coord_list(self.keywords['structure'].frac_coords,
                                       self.keywords['position'],
                                       atol=1e-04)

            struct_ed.replace_site(index, self.keywords['atom'])
        else:
            raise RuntimeError('Defect type %s not supported' % self.keywords['defecttype'])

        return struct_ed.modified_structure

    def _cart2frac(self, position):
        """Converts between cartesian coordinates and fractional coordinates"""
        return self.keywords['structure'].lattice.get_fractional_coords(self.keywords['position'])

    def write_input(self):
        """Writes the defected geometry to a file"""

    def _make_directory(self, directory):
        if os.path.exists(directory):
            print 'Directory %s exists!' % directory
        else:
            os.path.makedirs(directory)

