"""Tests for Recipetemplateparser"""

from MAST.parsers.recipetemplateparser import RecipeTemplateParser

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import InputOptions
from MAST.parsers import InputParser
from MAST.utility import MASTFile
import shutil

testname="recipeparser_test"
testdir = dirutil.get_test_dir(testname)

class TestRecipetemplateparser(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir(os.path.join(testdir,'workdir')):
            os.mkdir(os.path.join(testdir,'workdir'))

    def tearDown(self):
        if os.path.isdir(os.path.join(testdir,'workdir')):
            shutil.rmtree(os.path.join(testdir,'workdir'))

    def test___init__(self):
        raise SkipTest
        rt=os.path.join(testdir,'new_template.txt')
        wd=os.path.join(testdir,'workdir')
        iopt=InputOptions()
        pr=os.path.join(testdir,'workdir','personal_recipe.txt')
        myrtp=RecipeTemplateParser(inputOptions=iopt,working_directory=wd,templateFile=rt,personalRecipe=pr)
        self.assertTrue(os.path.isfile, myrtp.metafile)
        #self.testclass.__init__(**kwargs)

    def test_parse(self):
        #ipfile='small_workflow_with_scaling.inp'
        ipfile='phonon_with_neb.inp'
        workdir = os.path.join(testdir,'workdir')
        shutil.copy(os.path.join(testdir,ipfile),workdir)
        input_file = os.path.join(workdir,ipfile)
        parser_obj = InputParser(inputfile=input_file)
        input_options = parser_obj.parse()
        recipe_file_contents = input_options.get_item('recipe','recipe_file')
        myrtp = RecipeTemplateParser(inputOptions=input_options,
            working_directory=workdir,
            templateFile=recipe_file_contents,
            personalRecipe=input_file)
        parsed_lines = myrtp.parse()
        output_parsed = MASTFile()
        output_parsed.data = list(parsed_lines)
        output_parsed.to_file('workdir/personal_recipe_test_output.txt')
        self.assertEquals(len(myrtp.chunks),5)
        #mypr = MASTFile(pr)
        #for line in mypr.data:
        #    print line.rstrip()
        compare_pr = MASTFile(os.path.join(testdir,'compare','personal_recipe_test_output.txt'))
        self.assertEquals(parsed_lines, compare_pr.data)
        #self.testclass.parse()

    def test_parse_chunk(self):
        raise SkipTest
        #self.testclass.parse_chunk(chunk)

    def test_old_parsing(self):
        raise SkipTest
        #self.testclass.old_parsing()

    def test_process_system_name(self):
        raise SkipTest
        #self.testclass.process_system_name(processing_lines, system_name)

    def test_process_hop_combinations(self):
        raise SkipTest
        #self.testclass.process_hop_combinations(processing_lines, d_neblines)

    def test_process_images(self):
        raise SkipTest
        #self.testclass.process_images(processing_lines, n_images)

    def test_process_defects(self):
        raise SkipTest
        #self.testclass.process_defects(processing_lines, n_defects, d_defects)

    def test_process_phononlines(self):
        raise SkipTest
        #self.testclass.process_phononlines(processing_lines)

    def test_make_metadata_entries(self):
        raise SkipTest
        #self.testclass.make_metadata_entries(processing_lines)

    def test_get_unique_ingredients(self):
        raise SkipTest
        #self.testclass.get_unique_ingredients()

