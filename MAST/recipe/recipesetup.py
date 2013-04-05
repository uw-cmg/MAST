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

from MAST.ingredients.inducedefect import InduceDefect
from MAST.ingredients.optimize import Optimize
#from MAST.ingredients.singlepoint import SinglePoint
#from MAST.ingredients.optimizedefect import optimizeDefect

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import InputOptions
from MAST.utility import MASTError

from recipe_plan import RecipePlan
from pymatgen.core.structure import Structure


INGREDIENTS_LIBRARY = {\
                         'inducedefect'       : InduceDefect,\
                         'optimize'           : Optimize,\
#                         'singlepoint'        : SinglePoint,\
#                         'optimizedefect'     : OptimizeDefect,\
                      }

ALLOWED_KEYS = {
                  'recipeFile'     : (str, 'sic.recipe', 'Personalized recipe file'),\
                  'inputOptions'   : (InputOptions, None, 'Input options parsed using input parser'),\
                  'structure'      : (Structure, None, 'Structure to be used to create the ingredient objects'),\
               } 

DATA_PATH = "~/test_dir/"

class RecipeSetup(MASTObj):
    """Parses the personalized recipe file, create the recipe plan and prepare the ingredients
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.recipe_file    = self.keywords['recipeFile']
        self.input_options  = self.keywords['inputOptions']
        self.structure      = self.keywords['structure']
        self.program        = self.input_options.get_item('mast', 'program')
        self.scratch_dir    = self.input_options.get_item('mast', 'scratch_directory')
        self.delimiter      = '::'

        """Special keywords used in the recipe templates
        """
        self.recipe_keyword      = "recipe"
        self.ingredient_keyword  = "ingredient"
        self.parent_keyword      = "parent"
        self.child_keyword       = "child"

    def parse_recipe(self):
        """Parses the input personalized recipe file and takes 
           necessary info to setup the jobs
        """
        f_ptr            = open(self.recipe_file, "r")
        recipe_name      = ""
        ingredients_info = dict()

        for line in f_ptr.readlines():
            line = line.strip()
            #validate the line
            if not line or line.startswith('#'):
                continue

            #split line into elts
            elts         = [elt.strip() for elt in line.split()]
            init_keyword = elts[0].lower()

            #recipe name
            if init_keyword == self.recipe_keyword:
                recipe_name = elts[1]
                continue

            #handle job keyword
            if init_keyword == self.ingredient_keyword:
                ingredients_info.setdefault(elts[1], [elts[2], {}])
                continue

            #handle parent keyword
            if init_keyword == self.parent_keyword:
                parent_objs = []
                is_parent   = True
                for elt in elts[1:]:
                    if elt.lower() == self.child_keyword:
                        is_parent = False
                        continue
                    if is_parent:
                        parent_objs.append(elt)
                    else:
                        child_info      = elt.split(self.delimiter)
                        for parent in parent_objs:
                            if parent not in ingredients_info:
                                MASTError(self.__class__.__name__, "Parent Ingredient %s used in relationship without definition !!!" % parent)
                            ingredients_info[parent][1][child_info[0]] = child_info[1:]

        f_ptr.close()
        return ingredients_info, recipe_name

    def create_ingredient(self, name, ingredient_type, child_dict):
        """Creates the ingredient based on the ingredient type
        """
        if ingredient_type not in INGREDIENTS_LIBRARY:
            MASTError(self.__class__.__name__, "Ingredient %s not found !!!" % ingredient_type)
        ingredient_name = os.path.join(self.scratch_dir, name)
        return INGREDIENTS_LIBRARY[ingredient_type](name=ingredient_name, structure= self.structure, \
                                                    program=self.program,\
                                                    program_keys=self.input_options.get_item('ingredients', ingredient_type), child_dict=child_dict)

    def create_recipe_plan(self, ingredients_info, recipe_name):
        """Creates a recipe object which has the ingredients and dependency information
        """
        recipe_obj = RecipePlan(recipe_name)
        for name, (ingredient_type, child_dict) in ingredients_info.iteritems():
             ingredient_obj = self.create_ingredient(name, ingredient_type, child_dict)
             recipe_obj.add_ingredient(name, ingredient_obj)
             for child in child_dict.iterkeys():
                 recipe_obj.add_parent(child, name)
        return recipe_obj

    def prepare_ingredients(self, recipe_plan):
        """Prepare the ingredients
        """
        for ingredient_name, ingredient_obj in recipe_plan.ingredient_iterator():
            ingredient_obj.write_directory()
             

    def start(self):
        """Starts the setup process, parse the recipe file
           Use the input options and recipe info to 
           create directories and classes required
        """
        ingredients_info, recipe_name = self.parse_recipe()
        print ingredients_info
        recipe_plan                   = self.create_recipe_plan(ingredients_info, recipe_name)
        self.prepare_ingredients(recipe_plan)
        return recipe_plan

