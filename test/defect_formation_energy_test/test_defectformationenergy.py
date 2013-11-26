"""Tests for Defectformationenergy"""

from MAST.utility.defect_formation_energy.defectformationenergy import DefectFormationEnergy

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
import shutil

testname="defect_formation_energy_test"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)
oldarchive = os.getenv("MAST_ARCHIVE")
      
class TestDefectformationenergy(unittest.TestCase):

    def setUp(self):
        os.environ['MAST_ARCHIVE'] = os.path.join(testdir,'archive')
        os.chdir(testdir)
        if not os.path.isdir("dferesults"):
            os.mkdir("dferesults")

    def tearDown(self):
        os.environ['MAST_ARCHIVE'] = oldarchive
        if os.path.isdir("dferesults"):
            shutil.rmtree("dferesults")


    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(directory=None, plot_threshold=0.01)

    def test__calculate_defect_formation_energies(self):
        recipepath = os.path.join(testdir, 'archive','GaAs_defects_AsGa_recipe_defects_20131125T220427')
        mydfe = DefectFormationEnergy(directory=recipepath)
        os.chdir("dferesults")
        mydfe._calculate_defect_formation_energies() 
        self.assertEqual(True,True)
        #self.testclass._calculate_defect_formation_energies()

    def test_get_total_energy(self):
        raise SkipTest
        #self.testclass.get_total_energy(directory)

    def test_get_fermi_energy(self):
        raise SkipTest
        #self.testclass.get_fermi_energy(directory)

    def test_get_structure(self):
        raise SkipTest
        #self.testclass.get_structure(directory)

    def test_get_potential_alignment(self):
        raise SkipTest
        #self.testclass.get_potential_alignment(perf_dir, def_dir)

    def test_get_defect_formation_energies(self):
        raise SkipTest
        #self.testclass.get_defect_formation_energies()

    def test_defect_formation_energies(self):
        raise SkipTest
        #self.testclass.defect_formation_energies()

    def test_dfe(self):
        raise SkipTest
        #self.testclass.dfe()

    def test_print_table(self):
        raise SkipTest
        #self.testclass.print_table()

