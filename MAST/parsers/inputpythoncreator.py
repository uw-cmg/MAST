############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
from __future__ import print_function
import os

import numpy as np
import pymatgen as pmg

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import MAST2Structure
from MAST.utility.mastfile import MASTFile

ALLOWED_KEYS = {\
                 'input_options'    : (InputOptions, None, 'Input file name'),\
               }

class InputPythonCreator(MASTObj):
    """Creates a runnable python file from a *.inp file
        Attributes:
            self.name = name for created python file
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        if not 'mydir' in self.keywords.keys():
            mydir = os.getcwd()
        else:
            mydir == self.keywords['mydir']
        self.name=""
        mylines = self.print_input_options()
        inputpy = MASTFile()
        inputpy.data = mylines
        if self.name == "":
            self.name='input.py'
        inputpy.to_file(os.path.join(mydir, self.name))

    def addnewlines(self, linestoadd):
        """Appends newline to each string in a list
            Args:
                linestoadd <list of str>: lines to add newline to.
        """
        newline = "\n"
        newlist=list()
        for myline in linestoadd:
            try:
                myline = myline + newline
            except:
                raise MASTError(self.__class__.__name__,"Error with " + str(myline))
            newlist.append(myline)
        return newlist

    def print_input_options(self):
        inputoptions = self.keywords['input_options']
        pln=list()
        pln.append("from MAST.utility import InputOptions")
        pln.append("import numpy as np")
        pln.append("#MAST INPUT OPTIONS")
        pln.append("inputoptions = InputOptions()")

        for sectionname,secvalue in inputoptions.options.iteritems():
            pln.append(" ")
            pln.append("###################################")
            pln.append("#Section: %3s" % sectionname)
            pln.append("####################################")
            pln.append(" ")
            pln.extend(inputoptions.print_python_section('inputoptions.options',
                            sectionname))
            pln.append("inputoptions.options['"+sectionname+"'] = " + sectionname)
        pln.extend(self.print_create_structure_section())
        pln.extend(self.print_buffet_section())
        return self.addnewlines(pln)

    def print_create_structure_section(self):
        """Print the structure creation section.
            This is a little complex and may not belong here.

            Returns:
                pln <list of str>: list of lines to print out
        """
        pln=list()
        iops = self.keywords['input_options']
        strposfile = iops.get_item('structure','posfile')
        if strposfile is None:
            pln.append("iopscoords=inputoptions.get_item('structure','coordinates')")
            pln.append("iopslatt=inputoptions.get_item('structure','lattice')")
            pln.append("iopsatoms=inputoptions.get_item('structure','atom_list')")
            pln.append("iopsctype=inputoptions.get_item('structure','coord_type')")
            pln.append("from MAST.utility import MAST2Structure")
            pln.append("structure = MAST2Structure(lattice=iopslatt,")
            pln.append("    coordinates=iopscoords, atom_list=iopsatoms,")
            pln.append("    coord_type=iopsctype)")
        elif ('poscar' in strposfile.lower()):
            pln.append("from pymatgen.io.vaspio import Poscar")
            pln.append("structure = Poscar.from_file('" +strposfile+ "').structure")
        elif ('cif' in strposfile.lower()):
            from pymatgen.io.cifio import CifParser
            structure = CifParser(strposfile).get_structures()[0]
            pln.append("from pymatgen.io.cifio import CifParser")
            pln.append("structure = CifParser('"+strposfile+"').get_structures()[0]")
        else:
            error = 'Cannot build structure from file %s' % strposfile
            raise MASTError(self.__class__.__name__, error)
        pln.append("inputoptions.set_item('structure','structure',structure)")
        return pln

    def print_buffet_section(self):
        """Print the buffet section lines. Looping is accomplished
            through an "independent_loop_key" option key.

            Returns:
                pln <list of str>: list of lines to append
        """
        pln=list()
        pln.append(" ")
        pln.append("###############################################")
        pln.append("#Buffet Command Section")
        pln.append("################################################")
        pln.append(" ")
        pln.append("from MAST.recipe.recipeinput import RecipeInput")
        pln.append("from MAST.buffet.buffetmanager import Buffet")
        inputoptions = self.keywords['input_options']
        recipe_base_name = os.path.basename(inputoptions.options['recipe']['recipe_file'])
        recipe_base_name = recipe_base_name.split('.')[0]
        import time
        timestamp = time.strftime('%Y%m%d%H%M%S')
        self.name = recipe_base_name + '_' + timestamp + '.py'
        
        if "independent_loop_key" in inputoptions.options.keys():
            #for run_idx in xrange(len(SYSTEM_NAME)):
            #    recipe_input = RecipeInput(recipe_name="independent_recipe%s" % run_idx)

            #   #add input parameters
            #    recipe_input.addMastInfo(PROGRAM, SYSTEM_NAME[run_idx])
            #    recipe_input.addStructureInfoFromCoords(COORD_TYPE, COORDINATES[run_idx], LATTICE)
            #    recipe_input.addGlobalIngredientsInfo(ING_GLOBAL_PARAM)
            #    recipe_input.addIngredientsInfo(OPTIMIZE_ING, OPTIMIZE_INFO)
            #    recipe_input.addRecipeInfo(RECIPE_FILE)

            #    #add recipe input to buffet and cook it
            #    buffet_obj = Buffet(name="independent_recipe_%s" % SYSTEM_NAME[run_idx])
            #    buffet_obj.addRecipe(recipe_input)
            #    buffet_obj.cook()
            pass
        
        else:
            #copy from mast.py
            pln.append("from MAST.mast import MAST")
            pln.append("from MAST.ingredients.ingredients_loader import IngredientsLoader")
            pln.append("mast_obj = MAST()")
            pln.append("mast_obj.input_options = inputoptions")
            pln.append("ing_loader = IngredientsLoader()")
            pln.append("ing_loader.load_ingredients()")
            pln.append("ingredients_dict = ing_loader.ingredients_dict")

            #OMIT mast_obj._parse_input() here, because this is done above.
            pln.append("mast_obj._parse_recipe()")

            pln.append("mast_obj.initialize_environment()")
            pln.append("recipe_plan_obj = mast_obj._initialize_ingredients(ingredients_dict)")
            pln.append("mast_obj.pickle_plan(recipe_plan_obj)")
        return pln

