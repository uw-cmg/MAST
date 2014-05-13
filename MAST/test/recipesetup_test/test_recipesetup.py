"""Tests for Recipesetup"""

from MAST.recipe.recipesetup import RecipeSetup

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
from MAST.utility import MASTFile
import shutil
testname="recipesetup_test"
testdir = dirutil.get_test_dir(testname)

class TestRecipeSetup(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir(os.path.join(testdir,'workdir')):
            os.mkdir('workdir')

    def tearDown(self):
        if os.path.isdir(os.path.join(testdir,'workdir')):
            shutil.rmtree('workdir')

    def test___init__(self):
        rfile=os.path.join(testdir,'personalized')
        wdir=os.path.join(testdir,'workdir')
        struc=None
        iopt=InputOptions()
        myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        self.assertEqual(type(myrs.metafile),MAST.utility.metadata.Metadata)
        #self.testclass.__init__(**kwargs)

    def test_get_my_ingredient_options(self):
        raise SkipTest
        rfile=os.path.join(testdir,'personalized')
        wdir=os.path.join(testdir,'workdir')
        struc=None
        iopt=MAST.parsers.inputparser.InputParser(inputfile='input.inp').parse()
        myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        print iopt
        self.assertEqual(type(myrs.metafile),MAST.utility.metadata.Metadata)
        #self.testclass.get_my_ingredient_options(name, ingredient_type)

    def test_get_method_from_ingredient_type(self):
        raise SkipTest
        #self.testclass.get_method_from_ingredient_type(ingredtype, methodtype="")

    def test_create_recipe_plan(self):
        raise SkipTest
        #self.testclass.create_recipe_plan()

    def test_update_top_meta_for_ingred(self):
        raise SkipTest
        #self.testclass.update_top_meta_for_ingred(myingred)

    def test_create_ingredient(self):
        raise SkipTest
        #self.testclass.create_ingredient(my_ingred_input_options)

    def test_start(self):
        rfile=os.path.join(testdir,'personalized')
        wdir=os.path.join(testdir,'workdir')
        struc=None
        iopt=MAST.parsers.inputparser.InputParser(inputfile='input.inp').parse()
        myrs=RecipeSetup(recipeFile=rfile, workingDirectory=wdir, inputOptions=iopt, structure=struc)
        myrs.start()
        compare_metadata=MASTFile('compare_metadata.txt')
        myrs_metadata=MASTFile(os.path.join(wdir,'metadata.txt'))
        self.assertEqual(myrs_metadata.data, compare_metadata.data)
        dirlist=MAST.utility.dirutil.walkdirs(wdir)
        compare_dirlist=list()
        compare_dirlist.append('perfect_opt1')
        compare_dirlist.append('perfect_opt2')
        compare_dirlist.append('perfect_stat')
        compare_dirlist.append('inducedefect_group1')
        compare_dirlist.append('inducedefect_group3')
        compare_dirlist.append('defect_group1_q=n1_opt1')
        compare_dirlist.append('defect_group1_q=n1_opt2')
        compare_dirlist.append('defect_group1_q=n1_stat')
        compare_dirlist.append('phonon_group1_q=n1_int')
        compare_dirlist.append('phonon_group1_q=n1_int_parse')
        compare_dirlist.append('defect_group3_q=n1_opt1')
        compare_dirlist.append('defect_group3_q=n1_opt2')
        compare_dirlist.append('defect_group3_q=n1_stat')
        compare_dirlist.append('neb_group1-group3_q=n1_opt1')
        compare_dirlist.append('neb_group1-group3_q=n1_opt2')
        compare_dirlist.append('neb_group1-group3_q=n1_stat')
        compare_dirlist.append('phonon_group1-group3_q=n1_host1')
        compare_dirlist.append('phonon_group1-group3_q=n1_host1_parse')
        compare_dirlist_withpath=list()
        for diritem in compare_dirlist:
            compare_dirlist_withpath.append(os.path.join(wdir,diritem))
        mdiff=self.maxDiff
        self.maxDiff=None
        self.assertItemsEqual(dirlist,compare_dirlist_withpath)
        self.maxDiff=mdiff
        #self.testclass.start()

