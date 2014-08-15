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
from MAST.utility.finite_size_scaling.GenSC import genLMNs
from MAST.utility.finite_size_scaling.GenSC import gensc
from MAST.utility import dirutil
testname ="fss_gensc_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = os.path.join(dirutil.get_mast_install_path(),'test',testname)


class TestGenSC(unittest.TestCase):
    """Test GenSC
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass
    def test_genLMNs(self):
        primordial_struct = pymatgen.io.vaspio.Poscar.from_file("POSCAR_perfect").structure
        LMNlist = [[1, 1, 2], [2, 2, 2], [1, 2, 2], [1, 1, 1]]
        self.assertEqual(LMNlist,genLMNs(primordial_struct,minDefDist=5,maxNumAtoms=600,numStructAsked=5))
        
