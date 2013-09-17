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
import shutil

from MAST.utility import MASTObj
from MAST.utility import MAST2Structure
from MAST.utility import MASTError
from MAST.utility.picklemanager import PickleManager
from MAST.utility.dirutil import *
from MAST.utility import InputOptions

#from MAST.ingredients.ingredients_loader import IngredientsLoader

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

            self.origin_dir <str>: Original directory for mast -i *.inp command
            self.timestamp <str>: Timestamp for mast -i *.inp command
            self.asctime <str>: ASCII timestamp
            self.working_directory <str>: Working directory created for new recipe from mast -i *.inp command
            self.sysname <str>: System name (elements)
    """

    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options = None
        #self.recipe_name = None
        #self.structure = None
        #self.unique_ingredients = None
        self.origin_dir=""
        self.timestamp=""
        self.asctime=""
        self.working_directory=""
        self.sysname=""

    def check_independent_loops(self):
        """Checks for independent loops. If no independent loops are found,
            parse and run the input file. Otherwise, parse and run the
            generated set of input files.
        """
        ipl_obj = IndepLoopInputParser(inputfile=self.keywords['inputfile'])
        loopfiles = ipl_obj.main()
        #print "TTM DEBUG: loopfiles: ", loopfiles
        if len(loopfiles) == 0:
            self.parse_and_run_input()
        else:
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
        
        
        #run the *.py input script
        #import subprocess
        #oppath=self.input_options.get_item('mast','input_stem') + 'output'
        #opfile = open(oppath, 'ab')
        #run_input_script = subprocess.Popen(['python ' + ipc_filename], 
        #        shell=True, stdout=opfile, stderr=opfile)
        #run_input_script.wait()
        #opfile.close()
        self.start_from_input_options()
        return None
   

    def start_from_input_options(self):
        """Start the recipe template parsing and ingredient creation
            once self.input_options has been set.
        """
        print 'in start_from_input_options'
        
        #parse the recipe template file and create a personal file
        self.timestamp = time.strftime('%Y%m%dT%H%M%S')
        self.asctime = time.asctime()
        self.set_sysname()
        self.origin_dir = os.path.dirname(self.keywords['inputfile'])
        if self.origin_dir == "":
            self.origin_dir = os.getcwd()
        self.set_working_directory()

        self.make_working_directory()
        shutil.copy(self.keywords['inputfile'], os.path.join(self.working_directory,'input.inp'))
        self.parse_recipe_template()

        #make recipe plan object
        setup_obj = RecipeSetup(recipeFile=os.path.join(self.working_directory,'personal_recipe.txt'), 
                inputOptions=self.input_options,
                structure=self.input_options.get_item('structure','structure'), 
                workingDirectory=self.working_directory
                )
        recipe_plan_obj = setup_obj.start()
        recipe_plan_obj.print_status()

        self.pickle_plan(recipe_plan_obj)
        #self.pickle_input_options() 
        
        #create the *.py input script
        ipc_obj = InputPythonCreator(input_options=self.input_options)
        ipc_filename = ipc_obj.write_script()

    def make_working_directory(self):
        """Initialize the directory information for writing the 
            recipe and ingredient folders.
        """
        dir_path = self.working_directory
        try:
            os.mkdir(dir_path)
            topmeta = Metadata(metafile='%s/metadata.txt' % dir_path)
            topmeta.write_data('directory_created', self.asctime)
            topmeta.write_data('system_name', self.sysname)
            topmeta.write_data('origin_dir', self.origin_dir)
            topmeta.write_data('working_directory', self.working_directory)
            topmeta.write_data('timestamp', self.timestamp)
        except:
            MASTError(self.__class__.__name__, "Cannot create working directory %s !!!" % dir_path)


    def pickle_plan(self, recipe_plan_obj):
        """Pickles the reciple plan object to the respective file
           in the scratch directory
        """

        pickle_file = os.path.join(self.working_directory, 'mast.pickle')
        pm = PickleManager(pickle_file)
        pm.save_variable(recipe_plan_obj) 
   
    def pickle_input_options(self):
        """Temporary solution to input_options not being saved to the pickle correctly"""

        pickle_file = os.path.join(self.input_options.get_item('mast', 'working_directory'), 'input_options.pickle')
        pm = PickleManager(pickle_file)
        pm.save_variable(self.input_options)

    def parse_recipe_template(self):
        """Parses the recipe template file."""

        recipe_file = self.input_options.get_item('recipe', 'recipe_file')

        parser_obj = RecipeTemplateParser(templateFile=recipe_file, 
            inputOptions=self.input_options,
            personalRecipe=os.path.join(self.working_directory,'personal_recipe.txt'),
            working_directory=self.working_directory
            )
        self.input_options.update_item('recipe','recipe_name', parser_obj.parse())

    def set_input_stem_and_timestamp(self, input_options):
        """Set the input stem and timestamp.
            Args:
                input_options <InputOptions>
        """
        return
        mastkeys = input_options.get_section_keys('mast')
        if ('input_stem' in mastkeys) and ('timestamp' in mastkeys):
            return
        timestamp = time.strftime('%Y%m%dT%H%M%S')
        tstamp = time.asctime()
        inp_file = self.keywords['inputfile']
        stem_dir = os.path.dirname(inp_file)
        if len(stem_dir) == 0:
            stem_dir = dirutil.get_mast_scratch_path()
        inp_name = os.path.basename(inp_file).split('.')[0]
        stem_name = os.path.join(stem_dir, inp_name + '_' + timestamp + '_')
        input_options.update_item('mast', 'input_stem', stem_name)
        input_options.update_item('mast', 'timestamp', tstamp)

    def set_sysname(self):
        """Set system name."""
        element_str = self.get_element_string()
        system_name = self.input_options.get_item("mast", "system_name", "sys")
        self.sysname = system_name + '_' + element_str

    def set_working_directory(self):
        """Get the system name and working directory.
        """
        recipename = os.path.basename(self.input_options.get_item('recipe','recipe_file')).split('.')[0]
        dir_name = "%s_%s_%s" % (self.sysname, recipename, self.timestamp)
        dir_path = os.path.join(self.input_options.get_item('mast', 'scratch_directory'), dir_name)
        self.working_directory = dir_path

    def get_element_string(self):
        """Get the element string from the structure.
        """
        mystruc = self.input_options.get_item('structure','structure')
        elemset = set(mystruc.species)
        elemstr=""
        elemok=1
        while elemok == 1:
            try:
                elempop = elemset.pop()
                elemstr = elemstr + elempop.symbol
            except KeyError:
                elemok = 0
        return elemstr

