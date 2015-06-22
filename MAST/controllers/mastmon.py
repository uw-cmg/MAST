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
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.utility import loggerutils
from MAST.utility import InputOptions
from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.utility import MASTFile

class MASTMon(object):
    """The MAST monitor runs on a submission node and checks
        the status of each recipe in the MAST_SCRATCH directory.
        Attributes:
            self.scratch <str>: MAST_SCRATCH
            self._ARCHIVE <str>: MAST_ARCHIVE
            self.logger <logging logger>
    """ 
    def __init__(self):

        self.scratch = dirutil.get_mast_scratch_path()
        self._ARCHIVE = dirutil.get_mast_archive_path()
        self.make_directories() 
        self.logger = loggerutils.get_mast_logger('mast_monitor')
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.logger.info("\nMAST monitor started at %s.\n" % time.asctime())
        self.logger.info("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
        self.run(1)

    def make_directories(self):
        """Attempt to make scratch and archive directories
            if they do not exist.
        """
        try:
            if not os.path.exists(self.scratch):
                os.makedirs(self.scratch)
            if not os.path.exists(self._ARCHIVE):
                os.makedirs(self._ARCHIVE)
        except:
            raise MASTError(self.__class__.__name__,
                    "Error making directory for MASTmon and completed recipes")

    def check_recipe_dir(self, fulldir, verbose):
        """Check a recipe directory.
            Args:
                fulldir <str>: full path of recipe directory
                verbose <int>: verbosity
        """
        shortdir = os.path.basename(fulldir) #only the recipe directory name
        if not os.path.exists(fulldir):
            raise MASTError(self.__class__.__name__, "No recipe directory at %s" % fulldir)
        if os.path.exists(os.path.join(fulldir, "MAST_SKIP")):
            self.logger.warning("Skipping recipe %s due to the presence of a MAST_SKIP file in the recipe directory." % shortdir)
            return
        if os.path.exists(os.path.join(fulldir, "MAST_ERROR")):
            self.logger.error("ATTENTION!: Skipping recipe %s due to the presence of a MAST_ERROR file in the recipe directory." % shortdir)
            return
        self.logger.info("--------------------------------")
        self.logger.info("Processing recipe %s" % shortdir)
        self.logger.info("--------------------------------")
        os.chdir(fulldir) #need to change directories in order to submit jobs?
        myipparser = InputParser(inputfile=os.path.join(fulldir, 'input.inp'))
        myinputoptions = myipparser.parse()
        input_options_keys = myinputoptions.get_sections()
        key = 'personal_recipe'
        if key in input_options_keys:
            self.logger.debug("Key - personal recipe was found")
        personal_recipe_contents = myinputoptions.get_item('personal_recipe', 'personal_recipe_list')
        rsetup = RecipeSetup(recipeFile=personal_recipe_contents,
                inputOptions=myinputoptions,
                structure=myinputoptions.get_item('structure','structure'),
                workingDirectory=fulldir)
        recipe_plan_obj = rsetup.start()
        recipe_plan_obj.get_statuses_from_file()
        try:
            recipe_plan_obj.check_recipe_status(verbose)
        except Exception:
            import sys,traceback
            #ex_type, ex, trbck = sys.exc_info()
            errortext = traceback.format_exc()
            #del trbck
            errorfile = open(os.path.join(fulldir, "MAST_ERROR"), "ab")
            errorfile.write("ERROR LOGGED %s\n" % time.asctime())
            errorfile.write("%s\n" % errortext)
            errorfile.close()
            self.logger.warning("ERROR in recipe %s. Check MAST_ERROR file in the %s directory." % (shortdir, fulldir))
            #raise MASTError(self.__class__.__name__,"Error in recipe %s as follows: %s %s %s" % (shortdir, ex_type, ex, errortext))
        os.chdir(self.scratch)
        if recipe_plan_obj.status == "C":
            shutil.move(fulldir, self._ARCHIVE)
            summarypath = "%s/%s/SUMMARY.txt" % (self._ARCHIVE, shortdir)
            if os.path.isfile(summarypath):
                self.logger.info("Recipe %s completed." % shortdir)
                self.logger.info("SUMMARY.txt below:")
                summarytext = MASTFile(summarypath)
                for myline in summarytext.data:
                    self.logger.info(myline.strip())
        self.logger.info("-----------------------------")
        self.logger.info("Recipe %s processed." % shortdir)
        self.logger.info("-----------------------------")


    def run(self, verbose=0):
        """Run the MAST monitor.
        """
        curdir = os.getcwd()
        try:
            os.chdir(self.scratch)    
        except:
            os.chdir(curdir)
            errorstr = "Could not change directories to MAST_SCRATCH at %s" % self.scratch
            raise MASTError(self.__class__.__name__, errorstr)
        
        #dirutil.lock_directory(self.scratch, 1) # Wait 5 seconds
        #Directory is now locked by mast initially, but gets
        #unlocked at the end of the mastmon run.
        
        recipe_dirs = dirutil.walkdirs(self.scratch,1,1)
        if verbose == 1:
            self.logger.info("================================")
            self.logger.info("Recipe directories:")
            for recipe_dir in recipe_dirs:
                self.logger.info(recipe_dir)
            self.logger.info("================================")

        for recipe_dir in recipe_dirs:
            self.check_recipe_dir(recipe_dir, verbose)
                
        dirutil.unlock_directory(self.scratch) #unlock directory
        os.chdir(curdir)
