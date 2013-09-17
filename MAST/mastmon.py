import os
import time
import shutil
from MAST.utility.picklemanager import PickleManager
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup


class MASTmon(object):
    """MASTmon is a daemon to run dagscheduler class.
        This finds newly submitted recipe (recipe) and manage them. \n
        Also completed recipe is moved in archive directory by MASTmon.
        For consistency, users may touch sesssions in the archive directory."""
    
    def __init__(self):

        self.scratch = dirutil.get_mast_scratch_path()
        self._ARCHIVE = dirutil.get_mast_archive_path()
        self.make_directories() 

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
        curdir = os.getcwd()
        try:
            os.chdir(self.scratch)    
        except:
            os.chdir(curdir)
            errorstr = "Error: Failed to move to MASTmon scratch directory at %s" % self.scratch
            raise MASTError(self.__class__.__name__, errorstr)
        
        dirutil.lock_directory(self.scratch, 1) # Wait 5 seconds

        if verbose == 1:
            print "MAST is in: ", os.getcwd()
        
        recipe_dirs = dirutil.walkdirs(self.scratch,1,1)
        if verbose == 1:
            print "Recipe directories:"
            for recipe_dir in recipe_dirs:
                print recipe_dir

            # add new recipes
        for recipe_dir in recipe_dirs:
            self.check_recipe_dir(recipe_dir, verbose)
                
        dirutil.unlock_directory(self.scratch) #unlock directory
        os.chdir(curdir)
