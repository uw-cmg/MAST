"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import WriteIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile

testname="chopingredient_test"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestWriteIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_no_setup(self):
        raise SkipTest
        #self.testclass.no_setup()

    def test_write_neb(self):
        raise SkipTest
        #self.testclass.write_neb()

    def test_get_parent_structures(self):
        raise SkipTest
        #self.testclass.get_parent_structures()

    def test_get_parent_image_structures(self):
        raise SkipTest
        #self.testclass.get_parent_image_structures()

    def test_place_parent_energy_files(self):
        raise SkipTest
        #self.testclass.place_parent_energy_files()

    def test_write_neb_subfolders(self):
        raise SkipTest
        #self.testclass.write_neb_subfolders()

    def test_write_singlerun(self):
        raise SkipTest
        #self.testclass.write_singlerun()

    def test_write_phonon_multiple(self):
        raise SkipTest
        #self.testclass.write_phonon_multiple()

    def test_write_phonon_single(self):
        raise SkipTest
        #self.testclass.write_phonon_single()

    def test_get_my_phonon_params(self):
        raise SkipTest
        #self.testclass.get_my_phonon_params()

