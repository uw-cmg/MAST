import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
import shutil
import pymatgen
from pymatgen.io.vasp import Poscar
from pymatgen.analysis.ewald import EwaldSummation
from pymatgen.core.structure import Structure
import numpy as np
from MAST.utility.finite_size_scaling.EneVsVm import CalcV_M
from MAST.utility import dirutil
testname ="fss_enevsvm_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
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
        #V_M = round(5.185227,5)
        V_M = round(3.1339,5)
        struct = Poscar.from_file("POSCAR").structure
        dummystruct = Structure(struct.lattice,["F-"],[[0,0,0]])
        V_M_Ewald = EwaldSummation(dummystruct)
        print "Recip:", V_M_Ewald._recip
        print "Real:", V_M_Ewald._real
        print "Point:", V_M_Ewald._point #differs from previous version due to pymatgen removing incorrect jellium term (pymatgen commit a79a351398f2d34846f5eb8f682977bc85bad0d9)
        self.assertEqual(V_M,round(CalcV_M(struct),5)) 
