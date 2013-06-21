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
import glob
import inspect

from MAST import ingredients

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import InputOptions
from MAST.utility import MASTError

from recipe_plan import RecipePlan
from pymatgen.core.structure import Structure

ALLOWED_KEYS = {
                  'recipeFile'     : (str, None, 'Personalized recipe file'),\
                  'inputOptions'   : (InputOptions, None, 'Input options parsed using input parser'),\
                  'structure'      : (Structure, None, 'Structure to be used to create the ingredient objects'),\
                  'ingredientsDict': (dict, dict(), 'Dictionary of ingredients'),\
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
        self.ingredients_dict = self.keywords['ingredientsDict']
        self.delimiter      = '::'

        """Special keywords used in the recipe templates
        """
        self.recipe_keyword      = "recipe"
        self.ingredient_keyword  = "ingredient"
        self.parent_keyword      = "parent"
        self.child_keyword       = "child"
        print 'Setting up the recipe based on %s' % (self.recipe_file)

    def parse_recipe(self):
        """Parses the input personalized recipe file and takes 
           necessary info to setup the jobs
        """

        f_ptr            = open(self.recipe_file, "r")
        recipe_name      = ""
        ingredients_info = dict()

        for line in f_ptr.readlines():
            #print 'In parse_recipe:', line
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
                        if child_info[0] not in ingredients_info:
                            error = "Child Ingredient %s not defined!" % child_info[0]
                            raise MASTError(self.__class__.__name__, error)
                        for parent in parent_objs:
                            if parent not in ingredients_info:
                                error = "Parent Ingredient %s used in relationship without definition!" % parent
                                raise MASTError(self.__class__.__name__, error)
                            ingredients_info[parent][1][child_info[0]] = child_info[1:]

        f_ptr.close()
        return ingredients_info, recipe_name

    def create_ingredient(self, name, ingredient_type, child_dict):
        """Creates the ingredient based on the ingredient type.
            Notes:

            GRJ 6/19/2013: If the ingredient has not be specified in the input
                           file we create it here, and append on the defaults.
                           In addition, the defaults are appended on here, rather
                           than in InputParser (since each ingredient is made here
                           this makes more sense).
        """
        if ingredient_type not in self.ingredients_dict:
            error = "Ingredient '%s' requested by input file but not found in the recipe!" % ingredient_type
            raise MASTError(self.__class__.__name__, error)

        print '\nInitializing ingredient %s of type %s' % (name, ingredient_type)
        global_defaults = self.input_options.get_item('ingredients', 'global')

        if (ingredient_type not in self.input_options.get_section_keys('ingredients')):
            print 'Ingredient type %s has not be specified in the input file.' % ingredient_type
            print 'Using defaults from ingredients_global.'
            self.input_options.set_item('ingredients', ingredient_type, global_defaults)
        else:
            print 'Copying over defaults from ingredients_global for ingredient %s.' % ingredient_type
            ing_opt = self.input_options.get_item('ingredients', ingredient_type)
            for glob_key, glob_value in global_defaults.items():
                if glob_key not in ing_opt:
                    ing_opt[glob_key] = glob_value
            self.input_options.set_item('ingredients', ingredient_type, ing_opt)

        self.program = self.input_options.get_item('mast', 'program')
        self.scratch_dir = self.input_options.get_item('mast', 'working_directory')

        ingredient_name = os.path.join(self.scratch_dir, name)
        #print "TTM DEBUG: ",ingredient_type,":",self.input_options.get_item('ingredients',ingredient_type)
        #TTM update ingredients dict to include info from the 
        #'neb', 'defects', and 'chemical_potentials' sections
        pkey_d = self.input_options.get_item('ingredients', ingredient_type).copy()

        #print 'GRJ DEBUG: %s ingredient options %s' % (ingredient_type, pkey_d)
        #print 'GRJ DEBUG: %s ingredient' % ingredient_name
        #print 'GRJ DEBUG: Global ingredient options', self.input_options.get_item('ingredients', 'global')
        #print 'GRJ DEBUG: %s', self.input_options.options.keys()

        if 'defects' in self.input_options.options.keys():
            if 'defect_' in name.lower():
                pkey_d.update(self.input_options.get_item('defects','defects'))
                if 'inducedefect' not in ingredient_type:
                    clabel = [label for label in ingredient_name.split('_') if 'q=' in label][0].split('=')[1]
                    if 'n' in clabel[0]:
                        sign = -1
                    else:
                        sign = 1
                    charge = sign * int(clabel[1:])
                    #print 'GRJ DEBUG: Charge =', charge
                    pkey_d['charge'] = charge
                    #print 'Defect found, pkey_d =', pkey_d

        if 'neb' in self.input_options.options.keys():
            if 'neb' in name.lower():
                pkey_d.update(self.input_options.options['neb'])

        if 'chemical_potentials' in self.input_options.options.keys():
            pkey_d.update(self.input_options.options['chemical_potentials'])

        print 'Final pkey_d =', pkey_d
        return self.ingredients_dict[ingredient_type](name=ingredient_name, 
                    structure= self.structure, \
                    program=self.program, \
                    program_keys=pkey_d, \
                    child_dict=child_dict)

    def create_recipe_plan(self, ingredients_info, recipe_name):
        """Creates a recipe object which has the ingredients and dependency information
        """
        import inspect
        print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])

        recipe_obj = RecipePlan(recipe_name)

        for name, (ingredient_type, child_dict) in ingredients_info.iteritems():
             ingredient_obj = self.create_ingredient(name, ingredient_type, child_dict)
             recipe_obj.add_ingredient(name, ingredient_obj)

             for child in child_dict.iterkeys():
                 recipe_obj.add_parent(child, name)

        return recipe_obj

    def prepare_ingredients(self, recipe_plan):
        """Prepare the ingredients --- called after create_ingredients
        """
        #import inspect
        #print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])

        for ingredient_name, ingredient_obj in recipe_plan.ingredient_iterator():
            ingredient_obj.write_directory()

    def start(self):
        """Starts the setup process, parse the recipe file
           Use the input options and recipe info to 
           create directories and classes required
        """
        if self.recipe_file is None:
            raise MASTError(self.__class__.__name__, "Recipe file not provided!")
            
        if not os.path.exists(self.recipe_file):
            raise MASTError(self.__class__.__name__, "Recipe file not Found!")

        if not self.input_options:
            raise MASTError(self.__class__.__name__, "Input Options not provided!")

        if not self.ingredients_dict:
            raise MASTError(self.__class__.__name__, "Empty Ingredients Dict!")

        ingredients_info, recipe_name = self.parse_recipe()
        #print 'DEBUG:, ingredients info =', 
        #for ingredient, value in ingredients_info.items():
        #    print ingredient, value
        recipe_plan = self.create_recipe_plan(ingredients_info, recipe_name)
        self.prepare_ingredients(recipe_plan)
        return recipe_plan

