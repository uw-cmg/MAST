##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-11-11
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
                       'recdir' : (str, "", "Recipe directory"),\
                        'inputfile' : (str, "input.inp", "Input file")\
                   } 


class ModifyRecipe(MASTObj):
    """
        Modify a recipe according to the $recipe section of the
        input file:
            modify the $personal_recipe section
            add any folders that need to be added, with metadata, and status

        Attributes:
            self.recdir <str>: Recipe directory
            self.input_options <InputOptions object>: Stores the
                options parsed from the input file
            self.recipe_plan <RecipePlan object>: Stores the
                recipe plan parsed from self.input_options
            self.timestamp <str>: ASCII timestamp in YYYYMMDD'T'HHmmss format
            self.asctime <str>: ASCII timestamp
            self.logger <logging logger>
    """
    def __init__(self, **kwargs):
        MASTObj.__init__(self, ALLOWED_KEYS, **kwargs)
        self.recdir = self.keywords['recdir']
        self.input_options = None
        self.timestamp=""
        self.asctime=""
        self.recipe_plan = None
        self.logger = loggerutils.initialize_short_logger(os.path.join(get_mast_control_path(),"mast.log"))

    def modify_recipe(self):
        self.asctime = time.asctime()
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.logger.info("\nMAST modifying recipe started at %s using input file %s in directory %s" % (self.asctime, self.keywords['inputfile'], self.recdir))
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.check_for_structure_index_files()
        self.set_up_recipe()

    def check_for_structure_index_files(self):
        """Check for structure_index_files directory.
            If it exists, archive it.
        """
        if os.path.exists("%s/structure_index_files" % self.recdir):
            shutil.move("%s/structure_index_files" % self.recdir, "%s/archive_sif_%s" % (self.recdir, self.asctime))
            self.logger.info("Structure index file directory detected. Moved to an archive folder.")
        return

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
        self.update_recipe_metadata()
        self.parse_recipe_template()
        self.create_recipe_plan()
        self.create_archive_files()

    def create_recipe_plan(self):
        """Create the recipe plan object, and print its status.
        """
        personal_recipe_contents = self.input_options.get_item('personal_recipe', 'personal_recipe_list')
        setup_obj = RecipeSetup(recipeFile=personal_recipe_contents, 
                inputOptions=self.input_options,
                structure=self.input_options.get_item('structure','structure'), 
                workingDirectory=self.recdir
                )
        self.recipe_plan = setup_obj.start()
        self.recipe_plan.print_status()
        return



    def update_recipe_metadata(self):
        """Update the recipe metadata file.
        """
        topmeta = Metadata(metafile="%s/metadata.txt" % self.recdir)
        topmeta.write_data('recipe_updated', self.asctime)
        return

    def create_archive_files(self):
        """Save off archive files.
            Returns:
                creates archive_input_options.txt
                creates archive_recipe_plan.txt
        """
        inputsave = MASTFile()
        inputsave.data = repr(self.input_options)
        inputsave.to_file(os.path.join(self.recdir, 'archive_input_options_%s.txt' % self.timestamp))

        recipesave = MASTFile()
        recipesave.data = repr(self.recipe_plan)
        recipesave.to_file(os.path.join(self.recdir, 'archive_recipe_plan_%s.txt' % self.timestamp))

        return

    def parse_recipe_template(self):
        """Parses the recipe template file."""
        ipkeys = self.input_options.get_sections()
        if 'personal_recipe' in ipkeys:
            raise MASTError(self.__class__.__name__, "Personal recipe section still exists in input.inp. Remove this section and try the command again.")
            return
        recipe_file_contents = self.input_options.get_item('recipe', 'recipe_file')

        parser_obj = RecipeTemplateParser(templateFile=recipe_file_contents, 
            inputOptions=self.input_options,
            personalRecipe=os.path.join(self.recdir,'input.inp'),
            working_directory=self.recdir
            )
        personal_recipe_list = parser_obj.parse()
        #print personal_recipe_list
        if not personal_recipe_list:
            self.logger.info("Within mast/parse_recipe_template: Personal Recipe List is empty. Check whether input.inp has a personal_recipe section!")
        else:
            self.input_options.set_item('personal_recipe','personal_recipe_list',personal_recipe_list)

    def set_class_attributes(self):
        """Set class attributes, other than input options
        """
        time.sleep(1)
        self.timestamp = time.strftime('%Y%m%dT%H%M%S')
        self.asctime = time.asctime()

