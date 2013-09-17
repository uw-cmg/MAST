import os
from MAST.utility.picklemanager import PickleManager
#from MAST.DAG.dagscheduler import DAGScheduler
from MAST.utility import MASTError
from MAST.utility import dirutil
from MAST.parsers.inputparser import InputParser
from MAST.recipe.recipesetup import RecipeSetup
import time
#from MAST.DAG.dagutil import *
abspath = os.path.abspath
import shutil


class MASTmon(object):
    """MASTmon is a daemon to run dagscheduler class.
        This finds newly submitted recipe (recipe) and manage them. \n
        Also completed recipe is moved in archive directory by MASTmon.
        For consistency, users may touch sesssions in the archive directory."""
    
    def __init__(self):
        self.registered_dir = set()

        self.home = dirutil.get_mast_scratch_path()
        self._ARCHIVE = dirutil.get_mast_archive_path()

        self.scheduler = DAGScheduler()
        self.version = 0.1
        
        try:
            if not os.path.exists(self.home):
                os.makedirs(self.home)
            if not os.path.exists(self._ARCHIVE):
                os.makedirs(self._ARCHIVE)
        except:
            raise MASTError(self.__class__.__name__,
                    "Error making directory for MASTmon and completed recipes")

    def add_recipes(self, new_recipe_dirs, verbose):
        """recipe_dirs is a set of recipes in MASTmon home directory"""
        for recipe_dir in  new_recipe_dirs:
            fulldir = os.path.join(self.home, recipe_dir)
            #print 'recipe_dir =', recipe_dir
            if not os.path.exists(fulldir):
                raise MASTError("mastmon, add_recipes", "No recipe_dir at %s" % recipe_dir)
            #if not os.path.isfile(os.path.join(fulldir,'mast.pickle')):
            #    print "Skipping directory %s because there is no pickle file." % recipe_dir
            #    continue
            os.chdir(recipe_dir)
            #self.move_extra_files(fulldir)
            #if not os.path.isfile(os.path.join(fulldir, 'mast.pickle')):
            #    raise MASTError("mastmon, add_recipes", "No pickle file at %s/%s" % (fulldir, 'mast.pickle'))
            myipparser = InputParser(inputfile=os.path.join(fulldir, 'input.inp'))
            myinputoptions = myipparser.parse()
            rsetup = RecipeSetup(recipeFile=os.path.join(fulldir,'personal_recipe.txt'),
                    inputOptions=myinputoptions,
                    structure=myinputoptions.get_item('structure','structure'),
                    workingDirectory=fulldir)
            recipe_plan_obj = rsetup.start()

            #pm = PickleManager()
            #mastobj = pm.load_variable(os.path.join(fulldir,'mast.pickle'))
            recipe_plan_obj.get_statuses_from_file()
            recipe_plan_obj.check_recipe_status(verbose)
            #depdict = mastobj.dependency_dict
            #ingredients = mastobj.ingredients
            #pm.save(mastobj, os.path.join(fulldir, 'mast.pickle'))          
            os.chdir(self.home)
            if recipe_plan_obj.status == "C":
                shutil.move(fulldir, self._ARCHIVE)

            #if self.scheduler is None:
            #    print 'step 1: create DAGScheduler object'
            #    self.scheduler = DAGScheduler()
            
            #try: 
            #    self.scheduler.addjobs(ingredients_dict=ingredients, dependency_dict=depdict, sname=recipe_dir)    
            #except:
            #    raise MASTError(self.__class__.__name__,
            #        "Error adding jobs to scheduler.")
                
                
        #self.registered_dir = self.registered_dir.union(new_recipe_dirs)

    def del_recipe(self, scheduler, sid):
        print 'Deleting recipe %i from the scheduler' % sid
        scheduler.del_recipe(sid)

    def _save(self):
        """Save current stauts of MASTmon such as registered_dir and scheduler"""
        var_dict = {}
        var_dict['registered_dir'] = self.registered_dir
        var_dict['scheduler'] = self.scheduler
        var_dict['version']  = self.version
        self.pm.save(var_dict,filename=self.pn_mastmon)
        
    def _load(self):
        """Load MASTmon's information pickle file"""

        if os.path.isfile(self.pn_mastmon):
            var_dict = self.pm.load_variable(self.pn_mastmon)
            if 'version' in var_dict and var_dict['version'] != self.version:
                errorstr = "Error: mastmon_info.pickle is version %.2f while mastmon version is %.2f" % (var_dict['version'],self.version)
                raise MASTError(self.__class__.__name__, errorstr)

            if 'registered_dir' in var_dict:
                self.registered_dir = var_dict['registered_dir']
                
            if 'scheduler' in var_dict:
                self.scheduler = var_dict['scheduler']
        
    def run(self, niter=None, verbose=0, stopcond=None, interval=None, remove=None):
        """Run Mastmon. First of all, this makes MASTmon go to mastmon home load dagscheduler pickle.
            In addition, there are couple of options to run mastmon. \n
            ex) mastmon.run()  # run MASTmon forever as a real daemon. By default interval is 10 sec. \n
            ex) mastmon.run(interval=30) # run MASTmon forever as a real daemon. By default interval is 30 sec. \n
            ex) mastmon.run(niter=1) # run MASTmon one iteration for crontab user. By default interval is 10 sec. \n
            ex) mastmon.run(niter=20,stopcond='NOSESSION') # run MASTmon for 20 iterations. \n
            And stop it all recipes are done.
        """
        # move to mastmon home
        curdir = os.getcwd()
        try:
            os.chdir(self.home)    
        except:
            os.chdir(curdir)
            errorstr = "Error: Failed to move to MASTmon home %s" % self.home
            raise MASTError(self.__class__.__name__, errorstr)
        
        dirutil.lock_directory(self.home, 1) # Wait 5 seconds

        if verbose == 1:
            print "MAST is in: ", os.getcwd()
        if interval is None:
            interval = SCHEDULING_INTERVAL
            
        #load dagscheduler pickle
        #self._load()
        iter = 0;
        while True:
            if niter is not None and iter >= niter:
                break
            
            iter = iter + 1
            # get directories from mast home
            #recipe_dirs = os.walk('.').next()[1]

            #new_recipe_dirs = set(recipe_dirs) - self.registered_dir
            new_recipe_dirs = os.walk('.').next()[1]
            if verbose == 1:
                print "Recipe directories: ", new_recipe_dirs

            # add new recipes
            self.add_recipes(new_recipe_dirs, verbose)

            # run it for n iterations or until all recipes are complete
            #csnames = self.scheduler.run(niter=1, verbose=verbose)
            #self.scheduler.show_recipe_table()
            #remove complete recipes

            #if remove is not None:
            #    #print 'GRJ DEBUG: Before:', self.scheduler
            #    self.del_recipe(self.scheduler, remove)
            #    #print 'GRJ DEBUG: After:', self.scheduler
            #    self._save()
            #    break

            #self.registered_dir = self.registered_dir - csnames

            # save scheduler object
            #self._save()

            if stopcond is not None:
                if stopcond.upper() == 'NOSESSION' and len(self.registered_dir) == 0:
                    break
                
            #time.sleep(interval) #TTM remove this sleep
        # move back to original directory
        dirutil.unlock_directory(self.home) #unlock directory
        os.chdir(curdir)
            
                         
    def move_extra_files(self, recipedir):
        """Move extra files like input.py, output, and personalized recipe
            into the recipe directory.
            Args:
                recipedir <str>: Recipe directory
        """
        #mypm = PickleManager(os.path.join(self.home, recipedir, 'input_options.pickle'))
        myipparser = InputParser(inputfile=os.path.join(recipedir, 'input.inp'))
        myinputoptions = myipparser.parse()
        #myinputoptions = mypm.load_variable()
        workdir = myinputoptions.get_item("mast","working_directory")
        inpstem = myinputoptions.get_item("mast","input_stem")
        inpstembase = os.path.basename(inpstem)
                                     
        listdir = os.listdir(self.home)
        listdir.sort()
        for mystr in listdir:
            if inpstembase in mystr:
                os.rename(os.path.join(self.home, mystr), 
                    os.path.join(workdir, mystr))
            else:
                pass
                                     
        return
