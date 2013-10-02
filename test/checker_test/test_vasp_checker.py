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
        for fname in ["POSCAR","XDATCAR","DYNMAT","OSZICAR"]:
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

    def test_write_displacement(self):
        myvc = VaspChecker(name="dynamics")
        dp = myvc.read_my_displacement_file("dynamics")
        myvc.write_my_displacement_file("childdir", dp)
        dw = myvc.read_my_displacement_file("childdir")
        self.assertEqual(dp,dw)

    def test_forward_energy_file(self):
        myvc = VaspChecker(name="energy")
        myvc.forward_energy_file(os.path.join(testdir, "childdir"))
        op = pymatgen.io.vaspio.Outcar("energy/OSZICAR")
        oc = pymatgen.io.vaspio.Outcar("childdir/OSZICAR")
        self.assertEqual(op.run_stats,oc.run_stats)
    def test_is_complete(self):
        vcc = VaspChecker(name="done")
        vcs = VaspChecker(name="started")
        self.assertFalse(vcs.is_complete())
        self.assertTrue(vcc.is_complete())

    def test_is_ready(self):
        vcnr1 = VaspChecker(name="notready1")
        vcnr2 = VaspChecker(name="notready2")
        vcnr3 = VaspChecker(name="notready3")
        vcnr4 = VaspChecker(name="notready4")
        vcnr5 = VaspChecker(name="notready5")
        vcr = VaspChecker(name="ready")
        self.assertFalse(vcnr1.is_ready_to_run())
        self.assertFalse(vcnr2.is_ready_to_run())
        self.assertFalse(vcnr3.is_ready_to_run())
        self.assertFalse(vcnr4.is_ready_to_run())
        self.assertFalse(vcnr5.is_ready_to_run())
        self.assertTrue(vcr.is_ready_to_run())

    def test_combine_dynamical_matrix(self):
        raise SkipTest
        myvc = VaspChecker(name="dynamics_split")
        myvc.combine_dynamical_matrix_files(myvc.keywords['name'])

