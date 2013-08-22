############################################################################
# MAterials Simulation Toolbox (MAST)
# Version: January 2013
# Programmers: Tam Mayeshiba, Tom Angsten, Glen Jenness, Hyunwoo Kim,
#              Kumaresh Visakan Murugan
# Created at the University of Wisconsin-Madison.
# Replace this section with appropriate license text before shipping.
# Add additional programmers and schools as necessary.
############################################################################
import os

from MAST import ingredients

from MAST.utility import InputOptions
from MAST.utility import MASTObj
from MAST.utility import MASTError
from MAST.utility import Metadata

from pymatgen.core.structure import Structure
from MAST.recipe import recipeutility as ru
from MAST.recipe.recipeplan import RecipePlan
from MAST.ingredients.baseingredient import BaseIngredient
ALLOWED_KEYS = {
                  'recipeFile'     : (str, None, 'Personalized recipe file'),\
                  'inputOptions'   : (InputOptions, None, 'Input options parsed using input parser'),\
                  'structure'      : (Structure, None, 'Structure to be used to create the ingredient objects'),\
                  'ingredientsDict': (dict, dict(), 'Dictionary of ingredients'),\
               }

DATA_PATH = "~/test_dir/"

class RecipeSetup(MASTObj):
    """Parses the personalized recipe file,
        creates the recipe plan object,
        and prepares the ingredient directories.
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.recipe_file    = self.keywords['recipeFile']
        self.input_options  = self.keywords['inputOptions']
        self.structure      = self.keywords['structure']
        self.work_dir    = ""
        self.delimiter      = '::'

        self.metafile = Metadata(metafile='%s/metadata.txt' % self.input_options.get_item('mast', 'working_directory'))

        print 'Setting up the recipe based on %s' % (self.recipe_file)

    def get_my_ingredient_options(self, name, ingredient_type):
        """Creates the ingredient based on the ingredient type.
            Notes:

            GRJ 6/19/2013: If the ingredient has not be specified in the input
                           file we create it here, and append on the defaults.
                           In addition, the defaults are appended on here, rather
                           than in InputParser (since each ingredient is made here
                           this makes more sense).
            TTM 8/19/13: Unspecified ingredient is "default"
        """
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

        #self.program = self.input_options.options['ingredients'][ingredient_type]['mast_program']
        self.work_dir = self.input_options.get_item('mast', 'working_directory')

        ingredient_name = os.path.join(self.work_dir, name)
        pkey_d = self.input_options.get_item('ingredients', ingredient_type).copy()


        if 'defects' in self.input_options.options.keys():
            if 'defect_' in name.lower():
                for key in self.input_options.get_item('defects','defects'):
                    print key, ingredient_name

                pkey_d.update(self.input_options.get_item('defects','defects'))
                if 'inducedefect' not in ingredient_type:
                    data = self.metafile.read_data(ingredient_name.split('/')[-1]).split(',')
                    for datum in data:
                        if 'charge' in datum:
                            charge = int(datum.split(':')[-1])
                            pkey_d['mast_charge'] = charge

        if 'neb' in self.input_options.options.keys():
            if 'neb' in name.lower():
                pkey_d.update(self.input_options.options['neb'])
        if 'phonon' in self.input_options.options.keys():
            pkey_d.update(self.input_options.options['phonon'])
        if 'chemical_potentials' in self.input_options.options.keys():
            pkey_d.update(self.input_options.options['chemical_potentials'])

        allopt=dict()
        allopt['name'] = ingredient_name
        #allopt['program'] = self.program
        allopt['structure'] = self.structure
        allopt['program_keys'] = pkey_d
        return allopt
    def get_method_from_ingredient_type(self, ingredtype, methodtype=""):
        """Get write, ready, run, complete, and update methods from
            an ingredient type like volrelax_to_singlerun
            Args:
                ingredtype <str>: ingredient type
                methodtype <str>: method type ("mast_write_method", etc.)
            Returns:
                methodname <str>: method name ("write_default", etc.)
        """
        mdict=self.input_options.get_item("ingredients", ingredtype)
        if methodtype in mdict.keys():
            return mdict[methodtype]
        else:
            mdict=self.input_options.get_item("ingredients","global")
            if methodtype in mdict.keys():
                return mdict[methodtype]
            else:
                return None


    def create_recipe_plan(self):
        """Creates a recipe object which has the ingredients and dependency information
        """
        import inspect
        print 'GRJ DEBUG: %s.%s' % (self.__class__.__name__, inspect.stack()[0][3])
        [how_to_update, parents_to_check, how_to_run, recipe_name] = ru.read_recipe(self.recipe_file)

        recipe_obj = RecipePlan(recipe_name)
        ingredlist = how_to_run.keys()
        for ingred in ingredlist:
            self.update_top_meta_for_ingred(ingred)
            ingredtype = how_to_run[ingred]
            recipe_obj.ingredients[ingred]="I" #set all to initialized
            recipe_obj.ingred_input_options[ingred] = self.get_my_ingredient_options(ingred, ingredtype)
            self.create_ingredient(recipe_obj.ingred_input_options[ingred])
            recipe_obj.write_methods[ingred] = self.get_method_from_ingredient_type(ingredtype, "mast_write_method")
            recipe_obj.ready_methods[ingred] = self.get_method_from_ingredient_type(ingredtype, "mast_ready_method")
            recipe_obj.run_methods[ingred] = self.get_method_from_ingredient_type(ingredtype, "mast_run_method")
            recipe_obj.complete_methods[ingred] = self.get_method_from_ingredient_type(ingredtype, "mast_complete_method")
            recipe_obj.update_methods[ingred] = dict()
            for ichild in how_to_update[ingred].keys():
                updingredtype = how_to_update[ingred][ichild]
                recipe_obj.update_methods[ingred][ichild] = self.get_method_from_ingredient_type(updingredtype, "mast_update_children_method")
            recipe_obj.parents_to_check = dict(parents_to_check)
        return recipe_obj

    def update_top_meta_for_ingred(self, myingred):
        """Update the top metafile for an ingredient.
            Args:
                myingred <str>: ingredient name
        """
        datalist=list()
        datalist.append("ingredient type: %s " % myingred)
        if 'defect_' in myingred:
            defectlabel = myingred.split('defect_')[1].split('_')[0]
            if defectlabel.isdigit():
                defectlabel = "defect_" + defectlabel
            datalist.append("defect_label: %s" % defectlabel)
        if 'q=' in myingred:
            chargelabel = myingred.split('q=')[1].split('_')[0]
            datalist.append("charge: %s" % chargelabel)
        if 'neb_' in myingred:
            neblabel = myingred.split('neb_')[1].split('_')[0]
            datalist.append("neb_label: %s" % neblabel)
        if 'phonon_' in myingred:
            phononlabel = myingred.split('phonon_')[1].split('_')[0]
            datalist.append("phonon_label: %s" % phononlabel)
        data=','.join(datalist)
        self.metafile.write_data(myingred, data)

    def create_ingredient(self, my_ingred_input_options):
        """Create the ingredient directory and metadata file.
            Args:
                my_ingred_input_options <Input Options>: single ingredient
                                                            input options
        """
        allowed_keys = {
            'name' : (str, str(), 'Name of optimization directory'),
            'program': (str, str(), 'DFT program, e.g. "vasp"'),
            'program_keys': (dict, dict(), 'Dictionary of program keywords'),
            'structure': (Structure, None, 'Pymatgen Structure object')
            }


        myingred = BaseIngredient(allowed_keys, 
                        name=my_ingred_input_options['name'],
                        structure=my_ingred_input_options['structure'],
                        program_keys=my_ingred_input_options['program_keys'])
        myingred.write_directory()
        return


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


        #print 'DEBUG:, ingredients info =',
        #for ingredient, value in ingredients_info.items():
        #    print ingredient, value
        recipe_plan = self.create_recipe_plan()
        return recipe_plan

