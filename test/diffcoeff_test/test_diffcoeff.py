import pymatgen
import os
import shutil
import sys
import unittest
import time
import filecmp
from filecmp import dircmp
import subprocess

from MAST.ingredients.optimize import Optimize
from MAST.utility import MASTError
from pymatgen.io.vaspio import Poscar
from MAST.utility.dirutil import *
from nose import SkipTest

class TestDiffCoeff(unittest.TestCase):
    def setUp(self):
        scripts = get_mast_install_path()
        os.chdir(os.path.join(scripts, 'test','diffcoeff_test'))

        #self.myOpt=Optimize(name="test_opt1", program="vasp", program_keys={'ibrion':2, 'mast_kpoints':[4,4,4,"M"], 'mast_xc':"PBE", 'mast_setmagmom':"5 1", 'mast_adjustnelect':"-1"})

    def tearDown(self):
        pass
        #shutil.rmtree('test_opt1')
        #try:
        #    shutil.rmtree('test_opt2')
        #except OSError:
        #    pass
    
    def test_onefreq_run_from_prompt(self):
        #raise SkipTest
        verbose="0"
        mastpath = get_mast_install_path()
        myp=subprocess.Popen([mastpath+"/MAST/utility/diffusioncoefficient.py", mastpath+"/test/diffcoeff_test/diffcoeff_singlevac", "73", "1273", "100", "1", "w0=vac1-vac2",verbose])
            #stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #myp.communicate()[0]
        myp.wait()


    def test_fivefreq_run_from_prompt(self):
        verbose="1"
        mastpath = get_mast_install_path()
        myp=subprocess.Popen([mastpath+"/MAST/utility/diffusioncoefficient.py", mastpath+"/test/diffcoeff_test/alcu", "73", "1273", "100", "5", "w0=ep10-ep91,w1=ep10-ep37,w2=ep10-ep9,w3=ep10-ep20,w4=ep20-ep10",verbose])
            #stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #myp.communicate()[0]
        myp.wait()
