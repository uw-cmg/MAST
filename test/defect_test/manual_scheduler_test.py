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

def mastclear():
    dirlist = listscratch()
    for mydir in dirlist:
        shutil.rmtree(mydir)

def listscratch():
    scratchdir = os.getenv("MAST_SCRATCH")
    dirlist=dirutil.walkdirs(scratchdir, 1, 1, "*smalltestfcc*")
    print "SCRATCH DIRECTORY: ",dirlist
    return dirlist

def myrun(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return p

#mastclear()
mastrootdir = dirutil.get_mast_install_path()
inputfile = os.path.join(mastrootdir,'test/defect_test/mast.inp')
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
mastobj = pm.load_variable('mast.pickle')
depdict = mastobj.dependency_dict
ingredients = mastobj.ingredients

njobs = len(ingredients)

from MAST.DAG.dagscheduler import DAGScheduler
from MAST.DAG.jobentry import JobEntry
from MAST.DAG.sessionentry import SessionEntry
from MAST.DAG.jobentry import JobEntry
from MAST.DAG.jobtable import JobTable
from MAST.DAG.sessiontable import SessionTable
from MAST.DAG.dagutil import *

import time
# seession id
# Basic operation in DAGScheduler
# Default running mode 'noqsub' for machines which have no qsub
print 'step 1: create DAGScheduler object'
scheduler = DAGScheduler()

print 'step 2: add jobs'
print 'DEBUG:', ingredients
print 'DEBUG:', depdict
scheduler.addjobs(ingredients_dict=ingredients, dependency_dict = depdict)
print 'start scheduling'

TEST_CHECK_INTERVAL=10
scheduler._run_mode='serial'
while scheduler.has_incomplete_session():
    scheduler.run_jobs()
    time.sleep(TEST_CHECK_INTERVAL)
    scheduler.update_job_status()
    scheduler.show_session_table()

