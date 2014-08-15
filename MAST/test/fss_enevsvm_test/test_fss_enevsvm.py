import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
import shutil
import pymatgen
import numpy as np
from MAST.utility.finite_size_scaling.EneVsVm import CalcV_M
from MAST.utility import dirutil
testname ="fss_enevsvm_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = os.path.join(dirutil.get_mast_install_path(),'test',testname)


class TestEneVsVm(unittest.TestCase):
    """Test EneVsVm
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass
    def test_genLMNs(self):
        struct = pymatgen.io.vaspio.Poscar.from_file("POSCAR").structure
        V_M = round(5.185227,5)
        self.assertEqual(V_M,round(CalcV_M(struct),5)) 
