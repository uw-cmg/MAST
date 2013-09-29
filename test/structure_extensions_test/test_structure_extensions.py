import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
from MAST.ingredients.pmgextend.structure_extensions import StructureExtensions
import shutil
import pymatgen
import numpy as np

testname ="structure_extensions_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)


class TestSE(unittest.TestCase):
    """Test StructureExtensions
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass
    def test_induce_defect_frac(self):
        perfect = pymatgen.io.vaspio.Poscar.from_file("POSCAR_perfect").structure
        compare_vac1 = pymatgen.io.vaspio.Poscar.from_file("POSCAR_vac1").structure
        compare_int1 = pymatgen.io.vaspio.Poscar.from_file("POSCAR_int1").structure
        compare_sub1 = pymatgen.io.vaspio.Poscar.from_file("POSCAR_sub1").structure
        coord_type='fractional'
        threshold=1.e-4
        vac1={'symbol':'O', 'type': 'vacancy', 'coordinates':  np.array([0.25, 0.75, 0.25])}
        int1={'symbol':'Ni', 'type': 'interstitial', 'coordinates': np.array([0.3, 0.3, 0.3])}
        sub1={'symbol':'Fe', 'type': 'substitution','coordinates':np.array([0.25, 0.25,0.75])}
        sxtend = StructureExtensions(struc_work1=perfect)
        struc_vac1 = sxtend.induce_defect(vac1, coord_type, threshold)
        struc_int1 = sxtend.induce_defect(int1, coord_type, threshold)
        struc_sub1 = sxtend.induce_defect(sub1, coord_type, threshold)
        self.assertEqual(struc_vac1,compare_vac1)
        self.assertEqual(struc_int1,compare_int1)
        self.assertEqual(struc_sub1,compare_sub1)
