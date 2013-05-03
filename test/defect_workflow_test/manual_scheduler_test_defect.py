from MAST.utility.picklemanager import PickleManager
from MAST.parsers.recipeparsers import RecipeParser
from MAST.recipe.recipesetup import RecipeSetup
from MAST.recipe.recipe_plan import RecipePlan
import sys
import subprocess
import os
import copy
from MAST.utility import dirutil
import shutil
import time

def mastclear():
    dirlist = listscratch()
    for mydir in dirlist:
        shutil.rmtree(mydir)

def listscratch():
    scratchdir = os.getenv("MAST_SCRATCH")
    time.sleep(1)
    dirlist=dirutil.walkdirs(scratchdir, 1, 1, "*defecttest*")
    print "SCRATCH DIRECTORY: ",dirlist
    return dirlist

def myrun(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return p

mastclear()
mastrootdir = dirutil.get_mast_install_path()
inputfile = os.path.join(mastrootdir,'test/neb_workflow_test/defect_test.inp')
argv=[]
binpath = os.path.join(mastrootdir,'bin/mast')
bindir= os.path.join(mastrootdir,'bin')
print 'bindir : '+bindir
print 'binpath : '+binpath

print 'input directory : ' + inputfile
argv.append('-i')
argv.append(inputfile)
#sys.path.append(bindir)

subprocess.call([sys.executable, binpath, argv[0],argv[1]])

pm = PickleManager()
dirlist = listscratch()
topdir = ""
for mydir in dirlist:
    if len(mydir) > 0:
        topdir = mydir
print "TOP DIRECTORY: ",topdir
os.chdir(topdir)
sleepct=0
while not os.path.isfile('mast.pickle') and sleepct < 10:
    time.sleep(1)
    sleepct = sleepct + 1
if sleepct == 10:
    print "Timeout waiting for pickle."
    sys.exit()
mastobj = pm.load_variable('mast.pickle')
depdict = mastobj.dependency_dict
ingredients = mastobj.ingredients

njobs = len(ingredients)

from MAST.DAG.dagscheduler import DAGScheduler
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import DAGParser
from MAST.DAG.dagscheduler import SessionEntry
from MAST.DAG.dagscheduler import JobEntry
from MAST.DAG.dagscheduler import JobTable
from MAST.DAG.dagscheduler import SessionTable
from MAST.DAG.dagscheduler import set2str
from MAST.DAG.dagscheduler import JOB
import time
# seession id
# Basic operation in DAGScheduler
# Default running mode 'noqsub' for machines which have no qsub
print 'step 1: create DAGScheduler object'
scheduler = DAGScheduler()

print 'step 2: add jobs'
scheduler.addjobs(ingredients_dict=ingredients, dependency_dict = depdict)
print 'start scheduling'

TEST_CHECK_INTERVAL=10
scheduler._run_mode='serial'
while scheduler.has_incomplete_session():
    scheduler.run_jobs()
    time.sleep(TEST_CHECK_INTERVAL)
    scheduler.update_job_status()
    scheduler.show_session_table()

