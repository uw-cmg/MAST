#TTM from distribute_setup import use_setuptools
#TTM use_setuptools()
from setuptools.command.install import install
from setuptools import setup, find_packages
from distutils.command.build_py import build_py as _build_py

#python version check
import sys
import os
if sys.version_info[0] < 2 or (sys.version_info[0] == 2 and sys.version_info[1] < 7) or (sys.version_info[0] == 2 and sys.version_info[1] == 7 and sys.version_info[2] < 3):
    print "Python Version %d.%d.%d found" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    print "Python version >= 2.7.3 needed!!!"
    sys.exit(0)

#def run_initialize(my_platform):
#    import subprocess
#    myfile = open("setup_output","wb")
#    myinit = subprocess.Popen("mast_initialize", shell=True, stdout=myfile, stderr=myfile)
#    myinit.wait()
#    myfile.close()
#    return

def output_env_variable_info():
    myhome = os.getenv("HOME")
    print "==============================================="
    print "Add the following lines to your //home/user/.bashrc file"
    print "or a similar configuration file."
    print "Then, log out and log back in."
    print "See the MAST manual for more information."
    print "==============================================="

    print "export MAST_RECIPE_PATH=%s/MAST/recipe_templates" % myhome
    print "export MAST_SCRATCH=%s/MAST/SCRATCH" % myhome
    print "export MAST_ARCHIVE=%s/MAST/ARCHIVE" % myhome
    print "export MAST_CONTROL=%s/MAST/CONTROL" % myhome
    platform_list=""
    print "export MAST_PLATFORM=<choose one of: %s>" % platform_list
    return

def make_mast_tree(sourcedir):
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
        print "...Copying recipe directory from %s/recipe_templates into %s/MAST/recipe_templates" % (sourcedir,myhome)
        shutil.copytree("%s/recipe_templates" % sourcedir, "%s/MAST/recipe_templates" % myhome)
    else:
        print "...Directory %s/MAST/recipe_templates found; not creating." % myhome


class build_py(_build_py):
    """Specialized Python source builder."""
    print "Hello!"
    #make_mast_tree(_build_py.build_lib)
    output_env_variable_info()

setup(
        name="MAST_tam_test",
        packages=find_packages(),
        version="1.0.28",
        #setup_requires=["numpy>=1.6.1"],
        install_requires=["numpy>=1.6.1", "scipy>=0.10.1", "pymatgen>=2.8.8", "custodian>=0.5.1"],
        scripts=["MAST/bin/mast"],
        data_files=[
            ("MAST/recipe_templates",["files_to_copy/templates/neb_with_phonons.txt"])
        ],


        author="MAST Development Team, University of Wisconsin-Madison Computational Materials Group",
        author_email="ddmorgan@wisc.edu",
        #maintainer="Tam Mayeshiba",
        url="https://materialshub.org",
        license="MIT",
        description="MAterials Simulation Toolkit",
        long_description="MAterials Simulation Toolkit for diffusion and defects",
        keywords=["MAST","materials","simulation","diffusion","defects","ab initio","high throughput", "DFT", "density functional theory", "defect formation"],
        classifiers=[
            "Programming Language :: Python :: 2.7",
            "Development Status :: 4 - Beta",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
            "Topic :: Scientific/Engineering :: Information Analysis",
            "Topic :: Scientific/Engineering :: Physics",
            "Topic :: Scientific/Engineering :: Chemistry",
            "Topic :: Software Development :: Libraries :: Python Modules"
        ],
        cmdclass={'build_py': build_py}
)
