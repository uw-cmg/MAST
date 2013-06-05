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
from MAST.utility import MASTError
from MAST.parsers.recipeparsers import RecipeParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.recipe.recipeinput import RecipeInput

from MAST.ingredients.ingredients_loader import IngredientsLoader
from MAST.ingredients import *

ALLOWED_KEYS = {\
                   'recipe_input' : (RecipeInput, None, 'Recipe input instance'),\
               }


class RecipePlanner(MASTObj):
    """User interface to input an recipe input object and start planning
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs) 
        self.recipe_input = self.keywords['recipe_input']
        self.input_options = self.recipe_input.options

    def plan(self):
        """
            Creates the plan for the recipe and returns the plan object

            Args:
                None

            Returns:
                a recipe plan object which has the ingredients objects and the
                dependency dict between the ingredients

            Raises:
                None
        """ 
        #load ingredients
        ing_loader = IngredientsLoader()
        ing_loader.load_ingredients()
        ingredients_dict = ing_loader.ingredients_dict

        #parse the recipe
        self._parse_recipe()

        #initialize the environment and create the recipe plan
        self.initialize_environment()
        recipe_plan = self._initialize_ingredients(ingredients_dict)

        return recipe_plan


    def _parse_recipe(self):
        """
           Parses the generic recipe input file and gets the recipe name
           and the unique ingredients

           Args:
               None

           Returns:
               None

           Raises:
        """
        recipe_file = self.input_options.get_item('recipe', 'recipe_file')
        parser_obj = RecipeParser(templateFile=recipe_file, inputOptions=self.input_options, personalRecipe='test-recipe.txt')
        self.recipe_name = parser_obj.parse()
        self.unique_ingredients = parser_obj.get_unique_ingredients()


    def initialize_environment(self):
        """
           Create the working directory based on system name and timestamp
           and sets the input options

           Args:
               None

           Returns:
               None

           Raises:
        """
        system_name = self.input_options.get_item("mast", "system_name", "sys")
        #dir_name = "%s_%s_%s" % (system_name, self.recipe_name, time.strftime('%Y%m%dT%H%M%S'))
        dir_name = "%s_%s_%s_%s" % (system_name, self.recipe_name, self.input_options.get_item("main", "recipe_name"), time.strftime('%Y%m%dT%H%M%S'))
        dir_path = os.path.join(self.input_options.get_item('mast', 'scratch_directory'), dir_name)

        try:
            os.mkdir(dir_path)
        except:
            import pdb;pdb.set_trace()
            MASTError(self.__class__.__name__, "Cannot create working directory %s !!!" % dir_path)

        self.input_options.set_item('mast', 'working_directory', dir_path)


    def _initialize_ingredients(self, ingredients_dict):
        """
           Setup the recipe for the input options and create the recipe plan

           Args:
               ingredients_dict : a dict of ingredient names to the module object

           Returns:
               recipe plan instance created by the setup

           Raises:
               None
        """
        ingredient_global = self.input_options.get_item('ingredients', 'global')

        structure = self.input_options.get_item('structure', 'structure')

        for ingredient in self.unique_ingredients:
            if not self.input_options.get_item('ingredients', ingredient):
                self.input_options.set_item('ingredients', ingredient, ingredient_global)

        setup_obj = RecipeSetup(recipeFile='test-recipe.txt', inputOptions=self.input_options, structure=structure, ingredientsDict=ingredients_dict)
        recipe_plan = setup_obj.start()
        return recipe_plan

