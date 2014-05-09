##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-04-25
##############################################################
import os
import time
import shutil
import logging

from MAST.utility import MASTObj
from MAST.utility import MAST2Structure
from MAST.utility import MASTError
from MAST.utility import MASTFile
#from MAST.utility.picklemanager import PickleManager
from MAST.utility.dirutil import *
from MAST.utility import InputOptions
from MAST.utility import loggerutils
from MAST.parsers import InputParser
from MAST.parsers import IndepLoopInputParser
#from MAST.parsers import InputPythonCreator
from MAST.parsers.recipetemplateparser import RecipeTemplateParser
from MAST.recipe.recipesetup import RecipeSetup

from MAST.ingredients import *

ALLOWED_KEYS     = {\
                       'inputfile'    : (str, 'mast.inp', 'Input file name'),\
                       'outputfile'   : (str, 'mast.out', 'Output file name'),\
                   } 


class MAST(MASTObj):
    """User interface to set up a calculation group or a set of
        independently and/or pegged looped calculation groups.

        Attributes:
            self.input_options <InputOptions object>: Stores the
                options parsed from the input file
            self.recipe_plan <RecipePlan object>: Stores the
                recipe plan parsed from self.input_options and
                the personalized recipe file
            self.origin_dir <str>: Original directory for the
                                   mast -i *.inp command
            self.timestamp <str>: Timestamp of the mast -i *.inp
                                  command
            self.asctime <str>: ASCII timestamp
            self.working_directory <str>: Working directory 
                created for new recipe from mast -i *.inp command
            self.sysname <str>: System name (elements)
            self.logger <logging logger>
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.input_options = None
        self.origin_dir=""
        self.timestamp=""
        self.asctime=""
        self.working_directory=""
        self.sysname=""
        self.recipe_plan = None
        self.logger = loggerutils.initialize_short_logger(os.path.join(os.getenv("MAST_CONTROL"),"mast.log"))


    def check_independent_loops(self):
        """Checks for independent loops. If no independent loops are found,
            parse and run the input file. Otherwise, parse and run the
            generated set of input files.
        """
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.logger.info("\nMAST input started at %s with inputfile %s.\n" % (time.asctime(),self.keywords['inputfile']))
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        try:
            ipl_obj = IndepLoopInputParser(inputfile=self.keywords['inputfile'])
            loopfiles = ipl_obj.main()
            #print "TTM DEBUG: loopfiles: ", loopfiles
            if len(loopfiles) == 0:
                self.set_up_recipe()
            else:
                for ipfile in loopfiles:
                    self.keywords['inputfile']=ipfile
                    self.set_up_recipe()
        except MASTError as errormsg:
            self.logger.error(str(errormsg))
        self.logger.info("%%%%%%Finished processing inputfile.%%%%%%")


    def set_input_options(self):
        """Set input options.
        """
        parser_obj = InputParser(inputfile=self.keywords['inputfile'])
        self.input_options = parser_obj.parse()


    def set_up_recipe(self):
        """Set up the recipe.
            Set class attributes, 
            create the recipe directory,
            write recipe directory-level info,
            create the recipe plan (which creates ingredient
            directories), 
            and archive the recipe plan and input options.
        """
        self.set_input_options()
        self.set_class_attributes()
        self.make_working_directory()
        self.create_recipe_metadata()
        self.copy_input_file()
        self.parse_recipe_template()
        self.create_recipe_plan()
        self.create_archive_files()

    def create_recipe_plan(self):
        """Create the recipe plan object, and print its status.
        """
        setup_obj = RecipeSetup(recipeFile=os.path.join(self.working_directory,'personal_recipe.txt'), 
                inputOptions=self.input_options,
                structure=self.input_options.get_item('structure','structure'), 
                workingDirectory=self.working_directory
                )
        self.recipe_plan = setup_obj.start()
        self.recipe_plan.print_status()
        return

    def make_working_directory(self):
        """Initialize the directory for writing the 
            recipe and ingredient folders.
        """
        try:
            os.mkdir(self.working_directory)
        except:
            MASTError(self.__class__.__name__, "Cannot create working directory %s !!!" % self.working_directory)

    def copy_input_file(self):
        """Copy the input file to input.inp
            If the original file had loops, this is a single,
            non-looped input file.
        """
        shutil.copy(self.keywords['inputfile'], 
            os.path.join(self.working_directory,'input.inp'))

    def create_recipe_metadata(self):
        """Create the recipe metadata file.
        """
        topmeta = Metadata(metafile="%s/metadata.txt" % self.working_directory)
        topmeta.write_data('directory_created', self.asctime)
        topmeta.write_data('system_name', self.sysname)
        topmeta.write_data('origin_dir', self.origin_dir)
        topmeta.write_data('working_directory', self.working_directory)
        topmeta.write_data('timestamp', self.timestamp)
        return

    def create_archive_files(self):
        """Save off archive files.
            Returns:
                creates archive_input_options.txt
                creates archive_recipe_plan.txt
        """
        inputsave = MASTFile()
        inputsave.data = repr(self.input_options)
        inputsave.to_file(os.path.join(self.working_directory, 'archive_input_options.txt'))

        recipesave = MASTFile()
        recipesave.data = repr(self.recipe_plan)
        recipesave.to_file(os.path.join(self.working_directory, 'archive_recipe_plan.txt'))

        #pickle_plan = os.path.join(self.working_directory, 'archive_recipe_plan.pickle')
        #pm = PickleManager(pickle_plan)
        #pm.save_variable(self.recipe_plan)
        
        #pickle_options = os.path.join(self.working_directory, 'archive_input_options.pickle')
        #pm = PickleManager(pickle_options)
        #pm.save_variable(self.input_options)

        #create the *.py input script
        #ipc_obj = InputPythonCreator(input_options=self.input_options)
        #ipc_filename = ipc_obj.write_script(self.working_directory, 'archive_input_options.py')
        return

    def parse_recipe_template(self):
        """Parses the recipe template file."""

        recipe_file = self.input_options.get_item('recipe', 'recipe_file')

        parser_obj = RecipeTemplateParser(templateFile=recipe_file, 
            inputOptions=self.input_options,
            personalRecipe=os.path.join(self.working_directory,'personal_recipe.txt'),
            working_directory=self.working_directory
            )
        parser_obj.parse()

    def set_class_attributes(self):
        """Set class attributes, other than input options
        """
        time.sleep(1)
        self.timestamp = time.strftime('%Y%m%dT%H%M%S')
        self.asctime = time.asctime()
        self.set_sysname()
        self.set_origin_dir()
        self.set_working_directory()
    
    def set_origin_dir(self):
        """Set the origin directory (input file directory)
        """
        self.origin_dir = os.path.dirname(self.keywords['inputfile'])
        if self.origin_dir == "":
            self.origin_dir = os.getcwd()
        return

    def set_sysname(self):
        """Set system name."""
        element_str = self.get_element_map_string()
        if element_str == None:
            element_str = self.get_element_string()
        system_name = self.input_options.get_item("mast", "system_name", "sys")
        self.sysname = system_name + '_' + element_str
        return

    def set_working_directory(self):
        """Get the system name and working directory.
        """
        recipename = os.path.basename(self.input_options.get_item('recipe','recipe_file')).split('.')[0]
        #dir_name = "%s_%s_%s" % (self.sysname, recipename, self.timestamp)
        dir_name = "%s_%s" % (self.sysname, self.timestamp)
        dir_path = str(os.path.join(os.getenv("MAST_SCRATCH"), dir_name))
        self.working_directory = dir_path
        return

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
    def get_element_map_string(self):
        """Get the element string from the elementmap section.
            Returns:
                elemstr <str>: element string from elementmap
                                section of the input file,
                                in order of X1, X2, etc.
                If no elementmap section exists in the input file
                (which, when put in dictionary form, has the
                subsection name element_map), then the method
                returns None
        """
        elemmap = self.input_options.get_item('structure','element_map')
        if elemmap == None:
            return None
        elkeys = elemmap.keys()
        elkeys.sort()
        elstr=""
        for elkey in elkeys:
            elstr = elstr + elemmap[elkey]
        if elstr == "":
            return None
        return elstr


