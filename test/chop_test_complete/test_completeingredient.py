"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import IsCompleteIngredient

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile

testname="chop_test_complete"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestIsCompleteIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_complete_structure(self):
        raise SkipTest
        #self.testclass.complete_structure()

    def test_complete_singlerun(self):
        raise SkipTest
        #self.testclass.complete_singlerun()

    def test_complete_neb_subfolders(self):
        raise SkipTest
        #self.testclass.complete_neb_subfolders()

    def test_complete_subfolders(self):
        raise SkipTest
        #self.testclass.complete_subfolders()

