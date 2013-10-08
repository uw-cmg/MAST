"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import RunIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile

testname="chop_test_run"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestRunIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_run_singlerun(self):
        raise SkipTest
        #self.testclass.run_singlerun(mode='serial')

    def test_run_neb_subfolders(self):
        raise SkipTest
        #self.testclass.run_neb_subfolders()

    def test_run_subfolders(self):
        raise SkipTest
        #self.testclass.run_subfolders()

    def test_run_defect(self):
        raise SkipTest
        #self.testclass.run_defect()

    def test_run_strain(self):
        raise SkipTest
        #self.testclass.run_strain()

