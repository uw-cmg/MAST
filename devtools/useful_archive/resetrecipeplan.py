#!/usr/bin/env python
"""Reset the recipe plan pickle (mast.pickle) so that
    all ingredients show as "initialized"
    Args:
        [1] recipe path (can use ./ when inside recipe directory)
"""
from MAST.utility.picklemanager import PickleManager as pm
import sys
import os

if len(sys.argv) < 2:
    mypath = os.getcwd()
else:
    mypath = sys.argv[1]

fullpath = os.path.join(mypath, 'mast.pickle')
mypm = pm(fullpath)

my_recipe_plan = mypm.load_variable()
for ikey in my_recipe_plan.ingredients.keys():
    my_recipe_plan.ingredients[ikey] = "I"
mypm.save(my_recipe_plan,fullpath)
print "Reset all ingredient statuses to initialized."

