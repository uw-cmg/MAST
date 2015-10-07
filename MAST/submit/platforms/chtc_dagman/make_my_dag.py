##############################################################
# This code is part of the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2015-10-05
##############################################################
import sys
import os
import time
import shutil
import optparse # Allows for some command line option parsing
import glob
import subprocess 
import logging

from MAST.utility import MASTFile
from MAST.parsers import InputParser
from MAST.utility import dirutil
from MAST.utility import MASTError
from MAST.utility import loggerutils
from MAST.controllers.mastmon import MASTMon

def parse_arguments():
    if len(sys.argv) < 3:
        raise MASTError("DAG for MAST", "Not enough arguments.")
    recipe_name_raw = sys.argv[1]
    script_head_dir = sys.argv[2]
    recipe_name = recipe_name_raw
    sys.stdout.write("Parse arguments.\n")
    print "RECIPE: %s" % recipe_name
    #print_to_file(recipe_name, ing_name, "PARSED ARGUMENTS: %s" % (recipe_name, ing_name))
    return recipe_name, script_head_dir

def main():
    recipe_name, script_head_dir = parse_arguments()   
    mast_scratch = dirutil.get_mast_scratch_path()
    mymon = MASTMon()
    my_recipe_plan = mymon.set_up_recipe_plan(os.path.join(mast_scratch, recipe_name), 1)
    my_dag_contents=list()
    for iname in my_recipe_plan.ingredients: #all JOB lines need to be at top
        my_dag_contents.append("JOB %s submit.sh DIR %s\n" % (iname, iname))
    for iname in my_recipe_plan.ingredients:
        my_dag_contents.append("SCRIPT PRE %s %s/mast_do_setup.sh %s %s\n" % (iname, script_head_dir, recipe_name, iname)) 
        my_dag_contents.append("SCRIPT POST %s %s/mast_check_is_complete.sh %s %s\n" % (iname, script_head_dir, recipe_name, iname)) 
        my_dag_contents.append("RETRY %s 5\n" % iname)
        ptc = list(my_recipe_plan.parents_to_check[iname])
        for pname in ptc:
            my_dag_contents.append("PARENT %s CHILD %s\n" % (pname, iname))
    my_dag_file = MASTFile()
    my_dag_file.data = my_dag_contents
    my_dag_file.to_file(os.path.join(mast_scratch, recipe_name, "recipe.dag"))
    #print_to_file(recipe_name, ing_name, "MAIN: is complete: %s" % is_complete)
    return 0

if __name__ == '__main__':
    returncode = main()
    sys.exit(returncode)
