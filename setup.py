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


def output_env_variable_info():
    myhome = os.getenv("HOME")
    print ""
    print ""
    print ""
    print "==============================================="
    print "*************** ATTENTION *********************"
    print "Add the following lines to your //home/user/.bashrc file"
    print "or a similar configuration file."
    print "Then, log out and log back in."
    print "See the MAST manual for more information."
    print "==============================================="
    print ""
    print "export MAST_RECIPE_PATH=%s/MAST/recipe_templates" % myhome
    print "export MAST_SCRATCH=%s/MAST/SCRATCH" % myhome
    print "export MAST_ARCHIVE=%s/MAST/ARCHIVE" % myhome
    print "export MAST_CONTROL=%s/MAST/CONTROL" % myhome
    platform_list=["aci",
                    "bardeen",
                    "dlx",
                    "korczak",
                    "no_queue_system",
                    "pbs_generic",
                    "sge_generic",
                    "slurm_generic",
                    "stampede",
                    "turnbull"]
    myplatforms = "|".join(platform_list)
    print "export MAST_PLATFORM=<choose one of: %s>" % myplatforms
    print "==============================================="
    print ""
    print ""
    return

def make_mast_tree():
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


class build_py(_build_py):
    """Specialized Python source builder."""
    print "Starting setup for the MAterials Simulation Toolkit."
    make_mast_tree()
    output_env_variable_info()

setup(
        name="MAST_tam_test",
        packages=find_packages(),
        version="1.0.37",
        #setup_requires=["numpy>=1.6.1"],
        install_requires=["numpy>=1.6.1", "scipy>=0.10.1", "pymatgen>=2.8.8", "custodian>=0.5.1"],
        scripts=["MAST/bin/mast"],
        data_files=[
            ("MAST/recipe_templates",
                ["MAST/recipe_templates/neb_with_phonons.txt", 
                "MAST/recipe_templates/simple_optimization.txt",
                "MAST/recipe_templates/u_ramping.txt",
                "MAST/recipe_templates/defect_formation_energy.txt"]),
            ("MAST/examples",
                ["MAST/examples/README",
                "MAST/examples/neb_with_phonons.inp",
                "MAST/examples/simple_optimization.inp",
                "MAST/examples/u_ramping.inp",
                "MAST/examples/defect_formation_energy.inp",
                "MAST/examples/POSCAR.ga4as4"]),
            ("MAST/SCRATCH",
                ["MAST/initialization/README_Scratch"]),
            ("MAST/ARCHIVE",
                ["MAST/initialization/README_Archive"]),
            ("MAST/CONTROL",
                ["MAST/initialization/README_Control",
                "MAST/initialization/submitlist",
                "MAST/initialization/just_submitted"])
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
