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
    dirlist=dirutil.walkdirs(scratchdir, 1, 1, "*smalltestfcc*")
    print "SCRATCH DIRECTORY: ",dirlist
    return dirlist

def myrun(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    return p

#mastclear()
mastrootdir = dirutil.get_mast_install_path()
inputfile = os.path.join(mastrootdir,'test/neb_workflow_test/small_test.inp')
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

