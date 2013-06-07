############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan, Parker Sear
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os
import time

from MAST.utility import MASTObj
from MAST.utility import MAST2Structure
from MAST.utility import MASTError
from MAST.utility.picklemanager import PickleManager

from MAST.ingredients.ingredients_loader import IngredientsLoader

from MAST.parsers import InputParser
from MAST.parsers.recipeparsers import RecipeParser
from MAST.recipe.recipesetup import RecipeSetup

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
        self.recipe_name = None
        self.structure = None
        self.unique_ingredients = None
    
    def start(self):
        """Calls all the neccessary functions required to run MAST.
            Calls the following functions:
                _parse_input(): parses the various input files and 
                                fetches the options.
                _parse_recipe(): parses the recipe template file
        """

        ing_loader = IngredientsLoader()
        ing_loader.load_ingredients()
        ingredients_dict = ing_loader.ingredients_dict
        # print "Ingredients Dict : ", ingredients_dict

        self._parse_input()
        self._parse_recipe()

        self.initialize_environment()
        recipe_plan_obj = self._initialize_ingredients(ingredients_dict)
        self.pickle_plan(recipe_plan_obj)

    def initialize_environment(self):
        #mast_scratch_dir = os.cur_dir
        #if "MAST_SCRATCH_DIR" in os.environ:
        #    mast_scratch_dir = os.environ["MAST_SCRATCH_DIR"]
 
        system_name = self.input_options.get_item("mast", "system_name", "sys")

        dir_name = "%s_%s_%s" % (system_name, self.recipe_name, time.strftime('%Y%m%dT%H%M%S'))
        dir_path = os.path.join(self.input_options.get_item('mast', 'scratch_directory'), dir_name)

        try:
            os.mkdir(dir_path)
        except:
            MASTError(self.__class__.__name__, "Cannot create working directory %s !!!" % dir_path)

        self.input_options.set_item('mast', 'working_directory', dir_path)

    def pickle_plan(self, recipe_plan_obj):
        """Pickles the reciple plan object to the respective file
           in the scratch directory
        """

        pickle_file = os.path.join(self.input_options.get_item('mast', 'working_directory'), 'mast.pickle')
        pm = PickleManager(pickle_file)
        pm.save_variable(recipe_plan_obj) 
   
    def _parse_input(self):
        """Parses the input file"""
        parser_obj = InputParser(inputfile=self.keywords['inputfile'])
        self.input_options = parser_obj.parse()

    def _parse_recipe(self):
        """Parses the generic recipe file"""

        recipe_file = self.input_options.get_item('recipe', 'recipe_file')
        # print 'recipe_file =', recipe_file

        parser_obj = RecipeParser(templateFile=recipe_file, inputOptions=self.input_options,
                                  personalRecipe='test-recipe.txt')
        self.recipe_name = parser_obj.parse()

        self.unique_ingredients = parser_obj.get_unique_ingredients()

    def _initialize_ingredients(self, ingredients_dict):
        print '\nInitializing ingredients.'

        print '\nExtracting base structure.'
        structure = self.input_options.get_item('structure', 'structure')
        print structure, '\n'

        print '\nExtracting default ingredient options.'
        ingredient_global = self.input_options.get_item('ingredients', 'global')

        print '\nChecking status of defects.'
        ndefects = self.input_options.get_item('defects', 'num_defects')
        if ndefects == None:
            defects = False
            print 'No defects found.'
        else:
            defects = True
            defect_keys = self.input_options.get_item('defects', 'defects')
            if (ndefects == 1):
                print 'Found %i defect.\n' % ndefects
            else:
                print 'Found %i defects.\n' % ndefects

#        print "GRJ DEBUG:", self.unique_ingredients
#        print "GRJ DEBUG:", ingredients_dict
#        print "GRJ DEBUG:", defect_keys
        for ingredient in self.unique_ingredients:
            print 'Initializing ingredient %s.' % ingredient
            if (ingredient == 'inducedefect'):
                self.input_options.set_item('ingredients', ingredient, defect_keys)              
            elif (not self.input_options.get_item('ingredients', ingredient)):
                self.input_options.set_item('ingredients', ingredient, ingredient_global)

        print 'GRJ DEBUG:', self.input_options.get_item('ingredients', 'inducedefect')

        setup_obj = RecipeSetup(recipeFile='test-recipe.txt', inputOptions=self.input_options,
                                structure=structure, ingredientsDict=ingredients_dict)
        recipe_plan_obj = setup_obj.start()

        return recipe_plan_obj
