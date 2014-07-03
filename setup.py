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

myhome = os.getenv("HOME")
def output_env_variable_info():
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
    
class build_py(_build_py):
    """Specialized Python source builder."""
    print "Starting setup for the MAterials Simulation Toolkit."
    output_env_variable_info()

setup(
        name="MAST",
        packages=find_packages(),
        version="1.1.4",
        #setup_requires=["numpy>=1.6.1"],
        install_requires=["numpy>=1.6.1", "scipy>=0.10.1", "pymatgen>=2.8.8", "custodian>=0.5.1"],
        scripts=["MAST/bin/mast",
                "MAST/bin/mast_diffusion_coefficient",
                "MAST/bin/mast_defect_formation_energy"],
        data_files=[
            ("%s/MAST/recipe_templates" % myhome,
                ["MAST/recipe_templates/neb_with_phonons.txt", 
                "MAST/recipe_templates/simple_optimization.txt",
                "MAST/recipe_templates/u_ramping.txt",
                "MAST/recipe_templates/defect_formation_energy.txt"]),
            ("%s/MAST/examples" % myhome,
                ["MAST/examples/README",
                "MAST/examples/neb_with_phonons.inp",
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
            ("%s/MAST/CONTROL/platforms/aci" % myhome,
                ["MAST/submit/platforms/aci/mastmon_submit.sh",
                "MAST/submit/platforms/aci/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/bardeen" % myhome,
                ["MAST/submit/platforms/bardeen/mastmon_submit.sh",
                "MAST/submit/platforms/bardeen/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/dlx" % myhome,
                ["MAST/submit/platforms/dlx/mastmon_submit.sh",
                "MAST/submit/platforms/dlx/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/korczak" % myhome,
                ["MAST/submit/platforms/korczak/mastmon_submit.sh",
                "MAST/submit/platforms/korczak/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/no_queue_system" % myhome,
                ["MAST/submit/platforms/no_queue_system/mastmon_submit.sh",
                "MAST/submit/platforms/no_queue_system/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/pbs_generic" % myhome,
                ["MAST/submit/platforms/pbs_generic/mastmon_submit.sh",
                "MAST/submit/platforms/pbs_generic/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/sge_generic" % myhome,
                ["MAST/submit/platforms/sge_generic/mastmon_submit.sh",
                "MAST/submit/platforms/sge_generic/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/slurm_generic" % myhome,
                ["MAST/submit/platforms/slurm_generic/mastmon_submit.sh",
                "MAST/submit/platforms/slurm_generic/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/stampede" % myhome,
                ["MAST/submit/platforms/stampede/mastmon_submit.sh",
                "MAST/submit/platforms/stampede/submit_template.sh"]),
            ("%s/MAST/CONTROL/platforms/turnbull" % myhome,
                ["MAST/submit/platforms/turnbull/mastmon_submit.sh",
                "MAST/submit/platforms/turnbull/submit_template.sh"]),
            ("%s/MAST/CONTROL/citations" % myhome,
                ["MAST/summary/citations/mast_01",
                "MAST/summary/citations/mast_02",
                "MAST/summary/citations/pymatgen_01",
                "MAST/summary/citations/spglib_01",
                "MAST/summary/citations/structopt_01",
                "MAST/summary/citations/vasp_01",
                "MAST/summary/citations/vasp_02",
                "MAST/summary/citations/vasp_03",
                "MAST/summary/citations/vasp_04",
                "MAST/summary/citations/vaspneb_01",
                "MAST/summary/citations/vaspneb_02",
                "MAST/summary/citations/vaspneb_03",
                "MAST/summary/citations/vaspneb_04",
                "MAST/summary/citations/vaspneb_05",
                "MAST/summary/citations/vaspneb_06",
                "MAST/summary/citations/vasp_paw_01",
                "MAST/summary/citations/vasp_pps_01"]),
            ("%s/MAST/CONTROL/programkeys" % myhome,
                ["MAST/ingredients/programkeys/diff_allowed_keywords.py",
                "MAST/ingredients/programkeys/structopt_allowed_keywords.py",
                "MAST/ingredients/programkeys/vasp_allowed_keywords.py"]),
            ("%s/MAST/CONTROL" % myhome,
                ["MAST/structopt/Optimizer.py"])
        ],

        author="MAST Development Team, University of Wisconsin-Madison Computational Materials Group",
        author_email="ddmorgan@wisc.edu",
        #maintainer="Tam Mayeshiba",
        url="https://materialshub.org",
        license="MIT",
        description="MAterials Simulation Toolkit",
        long_description="MAterials Simulation Toolkit for diffusion and defects",
        keywords=["MAST","materials","simulation","diffusion","defects","ab initio","high throughput", "DFT", "density functional theory", "defect formation"],
        cmdclass={'build_py': build_py}
)
