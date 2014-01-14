#!/usr/bin/env python
import os
import subprocess
import shutil
print ""
print "Printing current MAST environment variables:"
print "MAST_INSTALL_PATH: ",os.getenv("MAST_INSTALL_PATH")
print "MAST_RECIPE_PATH: ",os.getenv("MAST_RECIPE_PATH")
print "MAST_SCRATCH: ",os.getenv("MAST_SCRATCH")
print "MAST_ARCHIVE: ", os.getenv("MAST_ARCHIVE")
print "MAST_CONTROL: ", os.getenv("MAST_CONTROL")
print "PYTHONPATH: ",os.getenv("PYTHONPATH")
print "PATH: ",os.getenv("PATH")
print "VASP_PSP_DIR: ",os.getenv("VASP_PSP_DIR")
print ""

print "Setting installation paths and new variable values:"
#Set $MAST_INSTALL_PATH, $PATH, and $PYTHONPATH
mycwd = os.getcwd()
mastinstall = "export MAST_INSTALL_PATH=%s" % mycwd
pypath = os.getenv("PYTHONPATH")
if pypath == None:
    exptpythonpath = "export PYTHONPATH=%s" % mycwd
else:
    exptpythonpath = "export PYTHONPATH=$PYTHONPATH:%s" % mycwd
mypath = os.getenv("PATH")
if mypath == None:
    exptpath = "export PATH=%s/bin" % mycwd
else:
    exptpath = "export PATH=$PATH:%s/bin" % mycwd

#Modify executables in $MAST_INSTALL_PATH/bin to be executable
modcmd = "chmod u+x %s/bin/*" % mycwd
modexec = subprocess.Popen(modcmd, shell=True)
modexec.wait()

#Make a MAST tree in user home
myhome = os.getenv("HOME")
print "...Making/looking for a MAST tree in %s/MAST..." % myhome
dirlist=list()
dirlist.append("%s/MAST" % myhome)
dirlist.append("%s/MAST/SCRATCH" % myhome)
dirlist.append("%s/MAST/ARCHIVE" % myhome)
dirlist.append("%s/MAST/CONTROL" % myhome)
for onedir in dirlist:
    if not os.path.isdir(onedir):
        print "...Creating directory %s" % onedir
        os.mkdir(onedir)
    else:
        print "...Directory %s found; not creating." % onedir

if not os.path.isdir("%s/MAST/recipe_templates" % myhome):
    print "...Copying recipe directory from %s/recipe_templates into %s/MAST/recipe_templates" % (mycwd,myhome)
    shutil.copytree("%s/recipe_templates" % mycwd, "%s/MAST/recipe_templates" % myhome)
else:
    print "...Directory %s/MAST/recipe_templates found; not creating." % myhome

filelist = list()
filelist.append("%s/MAST/CONTROL/submitlist" % myhome)
filelist.append("%s/MAST/CONTROL/just_submitted" % myhome)
for onefile in filelist:
    if not os.path.isfile(onefile):
        print "...Creating file %s" % onefile
        open(onefile, "ab").close()
    else:
        print "...File %s found; not creating." % onefile

#Get platform
print "Setting up platform."
platform_choices = list()
platform_choices.append("bardeen")
platform_choices.append("stampede")
platform_choices.append("dlx")
platform_choices.append("korczak")
platform_choices.append("slurm_generic")
platform_choices.append("pbs_generic")
plask = "What is your platform?\n\n"
plask = plask + "If your exact platform is not found (and even if it is), you may need to modify the files in $MAST_INSTALL_PATH/submit/platforms to conform to your platform.\n\n"
plask = plask + "Choose from:\n"
for platform in platform_choices:
    plask = plask + platform + '\n'
plask = plask + ">>> "

my_platform = "Not picked yet."
while not my_platform in platform_choices:
    my_platform = raw_input(plask)
    my_platform = my_platform.lower()
    if not my_platform in platform_choices:
        print "Platform choice %s not found. Please try again."
setfile = open("%s/MAST/CONTROL/set_platform" % myhome, "wb")
setfile.write(my_platform)
setfile.close()

shutil.copy("%s/submit/platforms/mastmon_submit_%s.sh" % (mycwd,my_platform), "%s/MAST/CONTROL/mastmon_submit.sh" % myhome)

#Print out environment variables
print "==============================================="
print "Add the following lines to your //home/user/.bashrc file"
print "or a similar configuration file."
print "Then, log out and log back in."
print "See the MAST manual for more information."
print "==============================================="

print mastinstall
print exptpythonpath
print exptpath
print "export MAST_RECIPE_PATH=%s/MAST/recipe_templates" % myhome
print "export MAST_SCRATCH=%s/MAST/SCRATCH" % myhome
print "export MAST_ARCHIVE=%s/MAST/ARCHIVE" % myhome
print "export MAST_CONTROL=%s/MAST/CONTROL" % myhome
