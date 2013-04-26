#!/usr/bin/env python
import os
mycwd = os.getcwd()
os.environ["MAST_INSTALL_PATH"]=mycwd
os.environ["MAST_RECIPE_PATH"]=os.path.join(mycwd, "recipe_templates")
vasppsp = os.getenv("VASP_PSP_DIR")
if vasppsp == None:
    print "Please set VASP_PSP_DIR to a pseudopotential directory."
print "MAST Environment Variables:"
print "MAST_INSTALL_PATH: ",os.getenv("MAST_INSTALL_PATH")
print "MAST_RECIPE_PATH: ",os.getenv("MAST_RECIPE_PATH")
print "VASP_PSP_DIR: ",os.getenv("VASP_PSP_DIR")
