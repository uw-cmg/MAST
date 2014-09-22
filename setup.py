########################################################################
# This is the setup script for the MAterials Simulation Toolkit (MAST)
# 
# Maintainer: Tam Mayeshiba
# Last updated: 2014-09-12
########################################################################
from setuptools.command.install import install
from setuptools import setup, find_packages
import sys
import os
import re


###Python version check
#print "Python version detected: %s" % sys.version_info
if sys.version_info[0] < 2 or (sys.version_info[0] == 2 and sys.version_info[1] < 7) or (sys.version_info[0] == 2 and sys.version_info[1] == 7 and sys.version_info[2] < 3):
    print "Python Version %d.%d.%d found" % (sys.version_info[0], sys.version_info[1], sys.version_info[2])
    print "Python version >= 2.7.1 needed!"
    sys.exit(0)

###Version load, adapted from http://stackoverflow.com/questions/2058802/how-can-i-get-the-version-defined-in-setup-py-setuptools-in-my-package/3619714#3619714
PKG = "MAST"
VERSIONFILE = os.path.join(PKG, "_version.py")
verstr = "unknown"
try:
    verstrline = open(VERSIONFILE, "rt").read()
except EnvironmentError:
    pass # Okay, there is no version file.
else:
    VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
    mo = re.search(VSRE, verstrline, re.M)
    if mo:
        verstr = mo.group(1)
    else:
        print "unable to find version in %s" % (VERSIONFILE,)
        raise RuntimeError("if %s.py exists, it is required to be well-formed" % (VERSIONFILE,))


###Set up home, and run setup

myhome = os.getenv("HOME")
setup(
        name="MAST",
        packages=find_packages(),
        version=verstr,
        install_requires=["numpy>=1.6.1", "scipy>=0.10.1", "pymatgen>=2.7.9", "custodian>=0.7.5"],
        scripts=["MAST/bin/mast",
                "MAST/bin/mast_diffusion_coefficient",
                "MAST/bin/mast_finite_size_scaling_sizes",
                "MAST/bin/mast_madelung_utility",
                "MAST/bin/mast_defect_formation_energy"],
        data_files=[
            ("%s/MAST/examples" % myhome,
                ["MAST/examples/README",
                "MAST/examples/neb_with_phonons.inp",
                "MAST/examples/finite_size_scaling.inp",
                "MAST/examples/simple_optimization.inp",
                "MAST/examples/u_ramping.inp",
                "MAST/examples/defect_formation_energy.inp",
                "MAST/examples/POSCAR.ga4as4"]),
            ("%s/MAST/SCRATCH" % myhome,
                ["MAST/initialization/README_Scratch"]),
            ("%s/MAST/ARCHIVE" % myhome,
                ["MAST/initialization/README_Archive"]),
            ("%s/MAST/CONTROL" % myhome,
                ["MAST/initialization/README_Control",
                "MAST/initialization/submitlist",
                "MAST/initialization/just_submitted",
                "MAST/submit/runmast.py"]),
        ],
        package_data={
            'MAST.submit.platforms.aci':['*.sh'],
            'MAST.submit.platforms.bardeen':['*.sh'],
            'MAST.submit.platforms.chtc':['*.sh','wrapper*','*mast*','test*','random*'],
            'MAST.submit.platforms.dlx':['*.sh'],
            'MAST.submit.platforms.korczak':['*.sh'],
            'MAST.submit.platforms.no_queue_system':['*.sh'],
            'MAST.submit.platforms.pbs_generic':['*.sh'],
            'MAST.submit.platforms.sge_generic':['*.sh'],
            'MAST.submit.platforms.slurm_generic':['*.sh'],
            'MAST.submit.platforms.stampede':['*.sh'],
            'MAST.submit.platforms.turnbull':['*.sh'],
            'MAST.utility.gbdiff.bin':['*'],
            'MAST.utility.gbdiff.data':['*'],
            'MAST.utility.gbdiff.doc':['*'],
            'MAST.utility.gbdiff.examples':['*'],
            'MAST.utility.gbdiff.middleware':['*'],
            'MAST.utility.gbdiff.rappture':['*'],
            'MAST.utility.gbdiff.src':['*'],
            'MAST.utility.diffanalyzer.bin':['*'],
            'MAST.utility.diffanalyzer.data':['*'],
            'MAST.utility.diffanalyzer.doc':['*'],
            'MAST.utility.diffanalyzer.examples':['*'],
            'MAST.utility.diffanalyzer.middleware':['*'],
            'MAST.utility.diffanalyzer.rappture':['*'],
            'MAST.utility.diffanalyzer.src':['*'],
            'MAST.summary.citations':['*'],
        },
        author="MAST Development Team, University of Wisconsin-Madison Computational Materials Group",
        author_email="ddmorgan@wisc.edu",
        #maintainer="Tam Mayeshiba",
        url="https://materialshub.org",
        license="MIT",
        description="MAterials Simulation Toolkit",
        long_description="MAterials Simulation Toolkit for diffusion and defects",
        keywords=["MAST","materials","simulation","diffusion","defects","ab initio","high throughput", "DFT", "density functional theory", "defect formation"],
        #cmdclass={'build_py': build_py}
)

print "***************************"
print "        ATTENTION!         "
print "***************************"
print "Please see the MAST documentation at http://pythonhosted.org/MAST"
print "in order to set up the following environment variables:"
print "MAST_SCRATCH"
print "MAST_ARCHIVE"
print "MAST_CONTROL"
print "MAST_PLATFORM"
print "(optionally, VASP_PSP_DIR)"

