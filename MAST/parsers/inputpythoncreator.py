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

ALLOWED_KEYS = {\
                 'input_options'    : (InputOptions, None, 'Input file name'),\
               }

class InputPythonCreator(MASTObj):
    """Creates a runnable python file from a *.inp file
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        if not 'mydir' in self.keywords.keys():
            mydir = os.getcwd()
        else:
            mydir == self.keywords['mydir']
        self.myfile = open(os.path.join(mydir, "input.py"), "wb")

    def print_input_options(self, lotsofspace=1):
        inputoptions = self.keywords['input_options']
        print("from MAST.utility import InputOptions", file=self.myfile)
        print("import numpy as np", file=self.myfile)
        print("#MAST INPUT OPTIONS", file=self.myfile)
        print("inputoptions = InputOptions()", file=self.myfile)

        for sectionname,secvalue in inputoptions.options.iteritems():
            if lotsofspace:
                print("", file=self.myfile)
                print("###################################", file=self.myfile)
            print("#Section: %3s" % sectionname, file=self.myfile)
            if lotsofspace:
                print("####################################", file=self.myfile)
                print("", file=self.myfile)
            #print "TTM DEBUG: straight print---------------------"
            #print secvalue
            #print "TTM DEBUG: New print:------------------------------"
            inputoptions.print_python_section('inputoptions.options',
                            sectionname, self.myfile)
        self.print_buffet_section()
        self.myfile.close()

    def print_buffet_section(self):
        """Print the buffet section lines. Looping is accomplished
            through an "independent_loop_key" option key.
        """
        print("", file=self.myfile)
        print("###############################################", 
                    file=self.myfile)
        print("#Buffet Command Section", file=self.myfile)
        print("################################################",
                    file=self.myfile)
        print("", file=self.myfile)
        print("from MAST.recipe.recipeinput import RecipeInput",
                    file=self.myfile)
        print("from MAST.buffet.buffetmanager import Buffet",
                    file=self.myfile)
        inputoptions = self.keywords['input_options']
        recipe_base_name = os.path.basename(inputoptions.options['recipe']['recipe_file'])
        recipe_base_name = recipe_base_name.split('.')[0]
        import time
        timestamp = time.strftime('%Y%m%dH%M%S')
        
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
            print("from MAST.mast import MAST",file=self.myfile)
            print("from MAST.ingredients.ingredients_loader import IngredientsLoader",file=self.myfile)
            print("mast_obj = MAST()",file=self.myfile)
            print("mast_obj.input_options = inputoptions",file=self.myfile)
            print("ing_loader = IngredientsLoader()",file=self.myfile)
            print("ing_loader.load_ingredients()",file=self.myfile)
            print("ingredients_dict = ing_loader.ingredients_dict",file=self.myfile)

            #OMIT mast_obj._parse_input() here, because this is done above.
            print("mast_obj._parse_recipe()",file=self.myfile)

            print("mast_obj.initialize_environment()",file=self.myfile)
            print("recipe_plan_obj = mast_obj._initialize_ingredients(ingredients_dict)",file=self.myfile)
            print("mast_obj.pickle_plan(recipe_plan_obj)",file=self.myfile)

