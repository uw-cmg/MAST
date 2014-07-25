"""Tests for Recipeutility"""

#from MAST.recipe.recipeutility import RecipeUtility

import unittest
from unittest import SkipTest
import os
import time
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTFile
from MAST.recipe import recipeutility
testname="recipeutility_test"
testdir = dirutil.get_test_dir(testname)

class TestRecipeUtility(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        pass

    def test_read_recipe(self):
        file_contents = MASTFile(os.path.join(testdir,'neb_int_personalized'))
        [htu, ptc,htr]=recipeutility.read_recipe(file_contents.data)
        #print htu
        #print ptc
        #print htr
        #print rname
        self.assertEqual(htu['neb_group1-group3_q=n1_opt2']['neb_group1-group3_q=n1_stat'],'neb_to_nebstat')
        self.assertItemsEqual(ptc['neb_group1-group3_q=n1_opt1'],['defect_group1_q=n1_stat','defect_group3_q=n1_stat'])
        self.assertEqual(htr['neb_group1-group3_q=n1_stat'],'nebstat_to_nebphonon')
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

