import os
import time
import shutil
from MAST.utility.picklemanager import PickleManager
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
import logging


class MASTmon(object):
    """MASTmon is a daemon to run dagscheduler class.
        This finds newly submitted recipe (recipe) and manage them. \n
        Also completed recipe is moved in archive directory by MASTmon.
        For consistency, users may touch sesssions in the archive directory."""
    
    def __init__(self):

        self.scratch = dirutil.get_mast_scratch_path()
        self._ARCHIVE = dirutil.get_mast_archive_path()
        self.make_directories() 
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        logger = logging.getLogger(__name__)
        logger.info("MAST monitor started at %s." % time.asctime())
        try:
            self.run(1)
        except BaseException as errormsg:
            logger.error(str(errormsg))

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
        if verbose > 0:
            print "Recipe directory: ", fulldir
        if not os.path.exists(fulldir):
            raise MASTError(self.__class__.__name__, "No recipe directory at %s" % fulldir)
        os.chdir(fulldir) #need to change directories in order to submit jobs?
        myipparser = InputParser(inputfile=os.path.join(fulldir, 'input.inp'))
        myinputoptions = myipparser.parse()
        rsetup = RecipeSetup(recipeFile=os.path.join(fulldir,'personal_recipe.txt'),
                inputOptions=myinputoptions,
                structure=myinputoptions.get_item('structure','structure'),
                workingDirectory=fulldir)
        recipe_plan_obj = rsetup.start()

        recipe_plan_obj.get_statuses_from_file()
        recipe_plan_obj.check_recipe_status(verbose)
        os.chdir(self.scratch)
        if recipe_plan_obj.status == "C":
            shutil.move(fulldir, self._ARCHIVE)


    def run(self, verbose=0):
        """Run the MAST monitor.
        """
        logging.basicConfig(filename="%s/mast.log" % os.getenv("MAST_CONTROL"), level=logging.DEBUG)
        logger=logging.getLogger(__name__)
        curdir = os.getcwd()
        try:
            os.chdir(self.scratch)    
        except:
            os.chdir(curdir)
            errorstr = "Could not change directories to MAST_SCRATCH at %s" % self.scratch
            raise MASTError(self.__class__.__name__, errorstr)
        
        dirutil.lock_directory(self.scratch, 1) # Wait 5 seconds
        
        recipe_dirs = dirutil.walkdirs(self.scratch,1,1)
        if verbose == 1:
            logger.info("Recipe directories:")
            for recipe_dir in recipe_dirs:
                logger.info(recipe_dir)
            logger.info("-------------------")

        for recipe_dir in recipe_dirs:
            logger.info("Processing recipe %s" % recipe_dir)
            try:
                self.check_recipe_dir(recipe_dir, verbose)
                logger.info("Recipe %s processed." % recipe_dir)
            except MASTError as errormsg:
                logger.error(str(errormsg))
                logger.info("Recipe %s failed with MASTError." % recipe_dir)
            except IOError as errormsg:
                logger.error(str(errormsg))
                logger.info("Recipe %s failed with IO Error." % recipe_dir)
                
        dirutil.unlock_directory(self.scratch) #unlock directory
        os.chdir(curdir)
