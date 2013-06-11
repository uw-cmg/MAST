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
from MAST.utility.dirutil import *

from MAST.ingredients.ingredients_loader import IngredientsLoader

from MAST.parsers import InputParser
from MAST.parsers import IndepLoopInputParser
from MAST.parsers import InputPythonCreator
from MAST.parsers.recipetemplateparser import RecipeTemplateParser
from MAST.recipe.recipesetup import RecipeSetup

from MAST.ingredients import *

ALLOWED_KEYS     = {\
                       'inputfile'    : (str, 'mast.inp', 'Input file name'),\
                       'outputfile'   : (str, 'mast.out', 'Output file name'),\
                   } 


class MAST(MASTObj):
    """User interface to set up a calculation group.

        Each instance of Interface sets up one calculation group.

        Attributes:
            self.input_options <InputOptions object>: used to store the options
                                           parsed from input file

    """

    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options = None
        #self.recipe_name = None
        #self.structure = None
        #self.unique_ingredients = None
   
    def check_independent_loops(self):
        """Checks for independent loops. If no independent loops are found,
            parse and run the input file. Otherwise, parse and run the
            generated set of input files.
        """
        ipl_obj = IndepLoopInputParser(inputfile=self.keywords['inputfile'])
        loopdict=ipl_obj.scan_for_indep_loop()
        if len(loopdict) == 0:
            self.parse_and_run_input()
        else:
            looplines=ipl_obj.prepare_looped_lines(loopdict)
            loopdata=ipl_obj.prepare_looped_datasets(looplines)
            loopfiles=ipl_obj.create_input_files(loopdata)
            for ipfile in loopfiles:
                self.keywords['inputfile']=ipfile
                self.parse_and_run_input()


    def parse_and_run_input(self):
        """
            Parses the *.inp input file and fetches the options.
            Sets an input stem name.
            Parses the recipe template file and creates a personalized recipe.
            Creates a *.py input script from the fetched options.
            Runs the *.py input script.
        """ 
        #parse the *.inp input file
        parser_obj = InputParser(inputfile=self.keywords['inputfile'])
        self.input_options = parser_obj.parse()
        
        #set an input stem name
        self.input_options.set_item('mast', 'input_stem', 
            self._initialize_input_stem())

        #parse the recipe template file and create a personal file
        self._parse_recipe_template()

        #create the *.py input script
        ipc_obj = InputPythonCreator(input_options=self.input_options)
        ipc_filename = ipc_obj.write_script()
        
        #run the *.py input script
        import subprocess
        run_input_script = subprocess.Popen(['python ' + ipc_filename], 
                shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE)
        run_input_script.wait()
        return None
   
    def _initialize_input_stem(self):
        inp_file = self.keywords['inputfile']
        if not ".inp" in inp_file:
            raise MASTError(self.__class__.__name__, "File '%s' should be named with a .inp extension and formatted correctly." % inp_file)
        
        timestamp = time.strftime('%Y%m%d%H%M%S')
        stem_dir = os.path.dirname(inp_file)
        if len(stem_dir) == 0:
            stem_dir = get_mast_scratch_path()
        inp_name = os.path.basename(inp_file).split('.')[0]
        stem_name = os.path.join(stem_dir, inp_name + '_' + timestamp + '_')
        return stem_name


    def start_from_input_options(self, input_options):
        """Start the recipe template parsing and ingredient creation
            once self.input_options has been set.
        """
        self.input_options = input_options
        ing_loader = IngredientsLoader()
        ing_loader.load_ingredients()
        ingredients_dict = ing_loader.ingredients_dict
        self.initialize_environment()
        recipe_plan_obj = self._initialize_ingredients(ingredients_dict)
        self.pickle_plan(recipe_plan_obj)


    def initialize_environment(self):
        #mast_scratch_dir = os.cur_dir
        #if "MAST_SCRATCH_DIR" in os.environ:
        #    mast_scratch_dir = os.environ["MAST_SCRATCH_DIR"]
 
        system_name = self.input_options.get_item("mast", "system_name", "sys")

        dir_name = "%s_%s_%s" % (system_name, self.input_options.get_item('recipe','recipe_name'), time.strftime('%Y%m%dT%H%M%S'))
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
   
    def _parse_recipe_template(self):
        """Parses the recipe template file."""

        recipe_file = self.input_options.get_item('recipe', 'recipe_file')

        parser_obj = RecipeTemplateParser(templateFile=recipe_file, 
                        inputOptions=self.input_options,
                        personalRecipe=self.input_options.get_item('mast','input_stem') + 'personal_recipe.txt')
        self.input_options.set_item('recipe','recipe_name', parser_obj.parse())
        self.unique_ingredients = parser_obj.get_unique_ingredients()

    def _initialize_ingredients(self, ingredients_dict):
        setup_obj = RecipeSetup(recipeFile=self.input_options.get_item('mast','input_stem') + 'personal_recipe.txt', 
                inputOptions=self.input_options,
                structure=self.input_options.get_item('structure','structure'), 
                ingredientsDict=ingredients_dict)
        recipe_plan_obj = setup_obj.start()
        return recipe_plan_obj
