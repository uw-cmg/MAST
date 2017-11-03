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
mysetuppath = os.path.dirname(os.path.abspath(__file__))
print "HOME: %s" % myhome
print "SETUP PATH: %s" % mysetuppath

def get_total_package_data():
    """Get total package data:
        * All files in MAST.submit.platforms.<platform name>
        * All files in MAST.summary.citation_files
        Returns:
            Dictionary with '<package name>':['<file name>','<file name>',...]
    """
    mydict = dict()
    #get platform files
    platformstem = os.path.join(mysetuppath, "MAST","submit","platforms")
    platformlist = os.listdir(platformstem)
    for platform in platformlist:
        myplatformdir = os.path.join(platformstem, platform)
        if os.path.isdir(myplatformdir):
            myfiles = os.listdir(myplatformdir)
            for myfile in myfiles:
                if os.path.isdir(os.path.join(myplatformdir, myfile)):
                    myfiles.remove(myfile)
            platformname = os.path.join("MAST","submit","platforms",platform).replace("/",".")
            mydict[platformname] = list(myfiles)
        else:
            pass
    #get citations files
    citationstem = os.path.join(mysetuppath,"MAST","summary","citation_files")
    citationlist = os.listdir(citationstem)
    mydict["MAST.summary.citation_files"] = list(citationlist)
    #get additional structopt files
    structoptstem = os.path.join(mysetuppath,"MAST","structopt")
    structoptlist=list()
    structoptlist.append("StructOpt_User_Guide_v1.docx")
    mydict["MAST.structopt"] = list(structoptlist)
    #get standalone files
    #mystandalones=['gbdiff','diffanalyzer']
    #for standalone in mystandalones:
    #    standstem = os.path.join(mysetuppath, "MAST","utility",standalone)
    #    standlist = os.listdir(standstem)
    #    mylist = list()
    #    for standitem in standlist:
    #        if os.path.isdir(os.path.join(standstem,standitem)):
    #            standfiles = os.listdir(os.path.join(standstem, standitem))
    #            for standitem2 in standfiles:
    #                mylist.append(os.path.join(standitem,standitem2))
    #        else:
    #            mylist.append(standitem)
    #    mydict[os.path.join("MAST","utility",standalone).replace("/",".")] = list(mylist)
    #get test files
    #teststem = os.path.join(mysetuppath, "MAST","test")
    #testlist = os.listdir(teststem)
    #for testname in testlist:
    #    mytestdir = os.path.join(teststem, testname)
    #    if os.path.isdir(mytestdir):
    #        mylist = list()
    #        myfiles = os.listdir(mytestdir)
    #        myfiles.sort()
    #        for myfile in myfiles:
    #            if os.path.isdir(os.path.join(mytestdir, myfile)):
    #                myfiles2 = os.listdir(os.path.join(mytestdir, myfile))
    #                myfiles2.sort()
    #                for myfile2 in myfiles2:
    #                    if os.path.isdir(os.path.join(mytestdir, myfile, myfile2)):
    #                        myfiles3 = os.listdir(os.path.join(mytestdir, myfile, myfile2))
    #                        myfiles3.sort()
    #                        for myfile3 in myfiles3:
    #                            if os.path.isdir(os.path.join(mytestdir, myfile, myfile2, myfile3)):
    #                                myfiles4 = os.listdir(os.path.join(mytestdir, myfile, myfile2, myfile3))
    #                                myfiles4.sort()
    #                                for myfile4 in myfiles4:
    #                                    mylist.append(os.path.join(myfile, myfile2, myfile3, myfile4))
    #                            else:
    #                                mylist.append(os.path.join(myfile, myfile2, myfile3))
    #                    else:
    #                        mylist.append(os.path.join(myfile, myfile2))
    #            else:
    #                mylist.append(myfile)
    #        mydict[os.path.join("MAST","test",testname).replace("/",".")]=mylist
    #    else:
    #        pass


    mykeys = mydict.keys()
    mykeys.sort()
    for key in mykeys:
        print key,": ", mydict[key]
    return mydict
setup(
        name="MAST",
        packages=find_packages(),
        version=verstr,
        install_requires=["numpy>=1.6.1", 
                "scipy>=0.10.1", 
                "pymatgen==4.7.1", 
                "custodian>=0.7.5",
                "pandas"],
        scripts=["MAST/bin/mast",
                "MAST/bin/mast_diffusion_coefficient",
                "MAST/bin/mast_finite_size_scaling_sizes",
                "MAST/bin/mast_madelung_utility",
                "MAST/utility/finite_size_scaling/mast_DFE_tool"],
        data_files=[
            ("%s/MAST/examples" % myhome,
                ["MAST/examples/README",
                "MAST/examples/neb_with_phonons.inp",
                "MAST/examples/neb_with_phonons_pathfinder.inp",
                "MAST/examples/neb_int_and_sub_test.inp",
                "MAST/examples/charged_defects_with_scaling.inp",
                "MAST/examples/simple_optimization.inp",
                "MAST/examples/u_ramping.inp",
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
        package_data=get_total_package_data(),
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

