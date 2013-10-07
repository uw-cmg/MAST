"""Tests for Chopingredient"""

from MAST.ingredients.chopingredient import UpdateChildrenIngredient

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

class TestUpdateChildrenIngredient(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test__fullpath_childname(self):
        raise SkipTest
        #self.testclass._fullpath_childname(childname)

    def test_give_structure(self):
        raise SkipTest
        #self.testclass.give_structure(childname)

    def test_give_neb_structures_to_neb(self):
        raise SkipTest
        #self.testclass.give_neb_structures_to_neb(childname)

    def test_give_saddle_structure(self):
        raise SkipTest
        #self.testclass.give_saddle_structure(childname)

    def test_give_phonon_multiple_forces_and_displacements(self):
        raise SkipTest
        #self.testclass.give_phonon_multiple_forces_and_displacements(childname)

    def test_give_phonon_single_forces_and_displacements(self):
        raise SkipTest
        #self.testclass.give_phonon_single_forces_and_displacements(childname)

    def test_give_structure_and_energy_to_neb(self):
        raise SkipTest
        #self.testclass.give_structure_and_energy_to_neb(childname)

    def test_give_structure_and_restart_files(self):
        raise SkipTest
        #self.testclass.give_structure_and_restart_files(childname)

