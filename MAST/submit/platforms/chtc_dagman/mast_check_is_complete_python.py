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

def print_to_file(recipe_name, ing_name, message):
    mypath = os.path.join(dirutil.get_mast_scratch_path(), recipe_name, ing_name, "DAG_post_script_output")
    with open(mypath, 'ab') as myfile:
        myfile.write(message)
        myfile.write("\n")
    return

def parse_arguments():
    if len(sys.argv) < 3:
        raise MASTError("DAG for MAST", "Not enough arguments.")
    recipe_name_raw = sys.argv[1]
    ing_name_raw = sys.argv[2]
    recipe_name = recipe_name_raw
    ing_name = ing_name_raw
    sys.stdout.write("Parse arguments.\n")
    print "RECIPE: %s" % recipe_name
    print "ING: %s" % ing_name
    #print_to_file(recipe_name, ing_name, "PARSED ARGUMENTS: %s" % (recipe_name, ing_name))
    return recipe_name, ing_name

def main():
    recipe_name, ing_name = parse_arguments()   
    mymon = MASTMon()
    my_recipe_plan = mymon.set_up_recipe_plan(recipe_name, 1)
    my_recipe_plan.ingred_to_check = ing_name
    is_complete = my_recipe_plan.complete_ingredient(ing_name)
    print "COMPLETE: %s" % is_complete
    #print_to_file(recipe_name, ing_name, "MAIN: is complete: %s" % is_complete)
    if is_complete:
        return 0
    return -1

if __name__ == '__main__':
    returncode = main()
    sys.exit(returncode)
