"""Tests for Recipe for FSS"""

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil

testname="fss_recipe_test"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestRecipeUtility(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test_read_recipe(self):
        [htu, ptc,htr, rname]=MAST.recipe.recipeutility.read_recipe(os.path.join(testdir,'fss_recipe.txt'))
        self.assertEqual(htu['defect_2x2x3_int1_opt']['neb_2x2x3_int1-int2_opt'],'relax_to_neb')
        self.assertItemsEqual(ptc['neb_1x2x4_int1-int2_opt'],['defect_1x2x4_int1_opt','defect_1x2x4_int2_opt'])
        self.assertEqual(htr['phonon_1x1x2_int1-int2_w0'],'phonon')
        #self.testclass.read_recipe(filename, verbose=0)

    def test_split_into_subrecipes(self):
        raise SkipTest
        #self.testclass.split_into_subrecipes(mydata)

    def test_get_indent_level(self):
        raise SkipTest
        #self.testclass.get_indent_level(myline)

    def test_parse_for_name_and_instructions(self):
        raise SkipTest
        #self.testclass.parse_for_name_and_instructions(myline)

    def test_parse_indentation_dict(self):
        raise SkipTest
        #self.testclass.parse_indentation_dict(idict)

    def test_make_indentation_dictionary(self):
        raise SkipTest
        #self.testclass.make_indentation_dictionary(subrlist)

