#!/usr/bin/env python
import os
mycwd = os.getcwd()
print "Current MAST Environment Variables:"
print "MAST_INSTALL_PATH: ",os.getenv("MAST_INSTALL_PATH")
print "MAST_RECIPE_PATH: ",os.getenv("MAST_RECIPE_PATH")
print "VASP_PSP_DIR: ",os.getenv("VASP_PSP_DIR")
print "MAST_SCRATCH: ",os.getenv("MAST_SCRATCH")
print "PYTHONPATH: ",os.getenv("PYTHONPATH")
print "PATH: ",os.getenv("PATH")

print "export MAST_INSTALL_PATH=" + mycwd
print "export MAST_RECIPE_PATH=" + mycwd + "/recipe_templates"
pypath = os.getenv("PYTHONPATH")
if pypath == None:
    print "export PYTHONPATH=" + mycwd
else:
    print "export PYTHONPATH=$PYTHONPATH:" + mycwd
mypath = os.getenv("PATH")
if mypath == None:
    print "export PATH=" + mycwd + "/bin"
else:
    print "export PATH=$PATH:" + mycwd + "/bin"
print "export MAST_SCRATCH=" + mycwd + "/SCRATCH"

print "env | grep MAST"
print "env | grep VASP"
