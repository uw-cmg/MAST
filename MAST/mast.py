############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from MAST.utility import MASTObj
from MAST.utility import MAST2Structure
from MAST.utility import MASTError

from MAST.parsers import InputParser
from MAST.parsers.recipeparsers import RecipeParser

from MAST.ingredients import *


#testing

ALLOWED_KEYS     = {\
                       'inputfile'    : (str, 'mast.inp', 'Input file name'),\
                       'outputfile'   : (str, 'mast.out', 'Output file name'),\
                   } 


class MAST(MASTObj):
    """User interface to set up a calculation group.

    Each instance of Interface sets up one calculation group.

    Attributes:
        options <InputOptions object>: used to store the options
                                       parsed from input file
    """

    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options = None
        self.structure = None
    
    def start(self):
        """Calls all the neccessary functions required to run MAST.
            Calls the following functions:
                _parse_input(): parses the various input files and fetches the options.
                _parse_recipe(): parses the recipe template file
        """
        self._parse_input()
        self._parse_recipe()
        self._initialize_ingredients()

    def _parse_input(self):
        """Parses the input file"""
        parser_obj = InputParser(inputfile=self.keywords['inputfile'])
        self.input_options = parser_obj.parse()

        self._build_structure()

        print self.input_options.get_item('defects', 'num_defects')

# Begin DEBUG section
#        from MAST.ingredients.inducedefect import InduceDefect as ID
#        for i in range(self.input_options.get_item('defects', 'num_defects')):
#            key = 'defect' + str(i)
#            print self.input_options.get_item('defects', 'defects')[key]
#
#            defect = self.input_options.get_item('defects', 'defects')[key]
#            print defect['coordinates'][None, :]
#            induced_defect = ID(atom=defect['symbol'],
#                                defecttype=defect['type'],
#                                position=tuple(defect['coordinates']),
#                                coordtype='fractional',
#                                structure=self.structure)
#            print induced_defect.induce_defect(), '\n'
# End DEBUG section

    def _build_structure(self):
        """Builds the structure object from input options
            For now we'll use the readers in pymatgen to generate the appropiate Structure object.

            For POSCAR and CIF files, code taken from the pymatget tutorial.
            Need to add flexibility for other file formats
        """
        posfile = self.input_options.get_item('structure', 'posfile')

        if (posfile is None):
            self.structure = MAST2Structure(lattice=self.input_options.get_item('unitcell', 'lattice'),
                                            coordinates=self.input_options.get_item('structure', 'coordinates'),
                                            atom_list=self.input_options.get_item('structure', 'atom_list'),
                                            coord_type=self.input_options.get_item('structure', 'coord_type'))
        elif ('poscar' in posfile.lower()):
            from pymatgen.io.vaspio import Poscar
            self.structure = Poscar.from_file(posfile).struct
        elif ('cif' in posfile.lower()):
            from pymatgen.io.cifio import CifParser
            self.structure = CifParser(posfile).get_structures()[0]
        else:
            error = 'Cannot build structure from file %s' % posfile
            MASTError(self.__class__.__name__, error)

# Begin DEBUG section
        print "\nPymatgen structure object:"
        print self.structure, '\n'
# End DEBUG section

    def _parse_recipe(self):
        """Parses the generic recipe file"""

        recipe_file = self.input_options.get_item('recipe', 'recipe_file')
        print 'recipe_file =', recipe_file

        parser_obj = RecipeParser(templateFile=recipe_file, inputOptions=self.input_options, personalRecipe='test-recipe.txt')
        parser_obj.parse()

    def _initialize_ingredients(self):
        print 'Initializing...'
        from MAST.ingredients import ingredient_dict as ID

