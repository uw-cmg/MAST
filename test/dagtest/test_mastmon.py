import os
import glob
import MAST.utility.picklemanger import PickleManger

cwd = os.getcwd()
pm = PickleManager()
pmastmon="mastmon_info.pickle"
completedir='archive'

Try:
    mastmon_home = os.environ['MAST_SCRATCH']
    os.chdir(mastmon_home)

    # Check if pickle of mastmon exists
    if os.path.isfile(pmastmon):
        mastmonobj = pm.load_variable(pmastmon)

    # Check only all top direcitories names
    recipe_directories = os.walk('.').next()[1]
    recipe_directories.remove(completedir)
    
    # Check mastmon sessions there is new directory which is not in session table

    # Go to that directory and add jobs

    # Come back to mastmon home
    
except:
    print 'Error'
    
os.chdir(Cwd)
