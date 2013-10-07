"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import IsReadyToRunIngredient

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

class TestIsReadyToRunIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_ready_singlerun(self):
        raise SkipTest
        #self.testclass.ready_singlerun()

    def test_ready_defect(self):
        raise SkipTest
        #self.testclass.ready_defect()

    def test_ready_neb_subfolders(self):
        raise SkipTest
        #self.testclass.ready_neb_subfolders()

    def test_ready_subfolders(self):
        raise SkipTest
        #self.testclass.ready_subfolders()

