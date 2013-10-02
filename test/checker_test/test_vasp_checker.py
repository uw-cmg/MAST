import os
import time
import unittest
from unittest import SkipTest
import filecmp
from filecmp import dircmp
import MAST
from MAST.ingredients.checker.vaspchecker import VaspChecker
import shutil
import pymatgen
import numpy as np

testname ="checker_test"
#oldcontrol = os.getenv("MAST_CONTROL")
#oldrecipe = os.getenv("MAST_RECIPE_PATH")
#oldscratch = os.getenv("MAST_SCRATCH")
#print "Old directories:"
#print oldcontrol
#print oldrecipe
#print oldscratch
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)


class TestVaspChecker(unittest.TestCase):
    """Test Vasp Checker
    """
    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        for fname in ["POSCAR","XDATCAR","DYNMAT"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass

    def test_get_structure_from_file(self):
        myvc = VaspChecker(name = "structure")
        getstr = myvc.get_structure_from_file("structure/POSCAR_perfect")
        mystr = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR_perfect").structure
        self.assertEqual(mystr, getstr)

    def test_get_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_initial_structure_from_directory()
        compare_nodir = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_initial_structure_from_directory("withdir")
        compare_withdir = pymatgen.io.vaspio.Poscar.from_file("withdir/POSCAR").structure
        self.assertEqual(withdir, compare_withdir)

    def test_get_final_structure(self):
        myvc = VaspChecker(name = "structure")
        nodir = myvc.get_final_structure_from_directory()
        compare_nodir = pymatgen.io.vaspio.Poscar.from_file("structure/CONTCAR").structure
        self.assertEqual(nodir, compare_nodir)
        withdir = myvc.get_final_structure_from_directory("withdir")
        compare_withdir = pymatgen.io.vaspio.Poscar.from_file("withdir/CONTCAR").structure
        self.assertEqual(withdir, compare_withdir)
    def test_forward_final_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_final_structure_file(os.path.join(testdir,"childdir"))
        fstrp = pymatgen.io.vaspio.Poscar.from_file("structure/CONTCAR").structure
        fstrc = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR").structure
        self.assertEqual(fstrp, fstrc)

    def test_forward_initial_structure(self):
        myvc = VaspChecker(name = "structure")
        myvc.forward_initial_structure_file(os.path.join(testdir,"childdir"))
        istrp = pymatgen.io.vaspio.Poscar.from_file("structure/POSCAR").structure
        istrc = pymatgen.io.vaspio.Poscar.from_file("childdir/POSCAR").structure
        self.assertEqual(istrp, istrc)
    
    def test_forward_dynamical_matrix(self):
        myvc = VaspChecker(name = "dynamics")
        myvc.forward_dynamical_matrix_file(os.path.join(testdir,"childdir"))
        dp = myvc.read_my_dynamical_matrix_file("dynamics")
        dc = myvc.read_my_dynamical_matrix_file("childdir")
        self.assertEqual(dp,dc)

    def test_forward_displacement(self):
        myvc = VaspChecker(name="dynamics")
        myvc.forward_displacement_file(os.path.join(testdir,"childdir"))
        dp = myvc.read_my_displacement_file("dynamics")
        dc = myvc.read_my_displacement_file("childdir")
        self.assertEqual(dp,dc)
