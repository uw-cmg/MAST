#TTM from distribute_setup import use_setuptools
#TTM use_setuptools()
from setuptools.command.install import install
from setuptools import setup, find_packages
from distutils.command.build_py import build_py as _build_py

#python version check
import sys
if sys.version_info[0] < 2 or (sys.version_info[0] == 2 and sys.version_info[1] < 7) or (sys.version_info[0] == 2 and sys.version_info[1] == 7 and sys.version_info[2] < 3):
    print "Python Version %d.%d.%d found" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    print "Python version >= 2.7.3 needed!!!"
    sys.exit(0)

def run_initialize(my_platform):
    import subprocess
    myfile = open("setup_output","wb")
    myinit = subprocess.Popen("mast_initialize", shell=True, stdout=myfile, stderr=myfile)
    myinit.wait()
    myfile.close()
    return

class build_py(_build_py):
    """Specialized Python source builder."""
    print "Hello!"
    my_platform = raw_input("Platform?: ")
    run_initialize(my_platform)

setup(
        name="MAST_tam_test",
        packages=find_packages(),
        version="1.0.13",
        #setup_requires=["numpy>=1.6.1"],
        install_requires=["numpy>=1.6.1", "scipy>=0.10.1", "pymatgen>=2.8.8", "custodian>=0.5.1"],
        scripts=["MAST/bin/mast","MAST/mast_initialize"],


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
        #data_files=[
        #    ("MAST/recipe_templates",["MAST/recipe_templates/neb_with_phonons.txt"])
        #],
        cmdclass={'build_py': build_py}
)
