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
        pass
        #if os.path.isdir(os.path.join(testdir,'workdir')):
        #    shutil.rmtree(os.path.join(testdir,'workdir'))

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
        #rt=os.path.join(testdir,'new_template.txt')
        rt=os.path.join(testdir,'phonon_with_neb.inp')
        wd=os.path.join(testdir,'workdir')
        parser_obj = InputParser(inputfile=rt)
        input_options = parser_obj.parse()
        iopt=InputOptions()
        iopt.options['neb']=dict()
        iopt.options['neb']['nebs']=dict()
        iopt.options['neb']['nebs']['group1-group2']=dict()
        iopt.options['neb']['nebs']['group1-group2']['phonon']=dict()
        iopt.options['neb']['nebs']['group1-group2']['phonon']['solvent']=dict()
        iopt.options['neb']['nebs']['group1-group2']['phonon']['solute']=dict()
        iopt.options['neb']['nebs']['group1-group3']=dict()
        iopt.options['neb']['nebs']['group1-group3']['phonon']=dict()
        iopt.options['defects']=dict()
        iopt.options['defects']['defects']=dict()
        iopt.options['defects']['defects']['group1']=dict()
        iopt.options['defects']['defects']['group1']['charge']=[-1,1,2,3]
        iopt.options['defects']['defects']['group1']['phonon']=dict()
        iopt.options['defects']['defects']['group2']=dict()
        iopt.options['defects']['defects']['group2']['charge']=[1,2,3]
        iopt.options['defects']['defects']['group2']['phonon']=dict()
        iopt.options['defects']['defects']['group3']=dict()
        iopt.options['defects']['defects']['group3']['charge']=[-1]
        iopt.options['defects']['defects']['group3']['phonon']=dict()
        iopt.options['defects']['defects']['group3']['phonon']['host1']=dict()
        recipe_file_contents = input_options.get_item('recipe', 'recipe_file')
        shutil.copy(rt,os.path.join(wd,'input.inp'))
        pr=os.path.join(testdir,'workdir','input.inp')
        myrtp=RecipeTemplateParser(inputOptions=iopt,working_directory=wd,templateFile=recipe_file_contents,personalRecipe=pr)
        rname=myrtp.parse()
        prf=os.path.join(testdir,'Personal_Recipe_File.txt')
        per_rec_file = open(prf, 'w')
        print "SASWATI DE *********************************"
        for line in rname:
            print line.rstrip()
            per_rec_file.write("%s\n" % line.rstrip())
        print "SASWATI DE *********************************"
        self.assertEquals(len(myrtp.chunks),5)
        #mypr = MASTFile(pr)
        #for line in mypr.data:
        #    print line.rstrip()
        compare_pr = MASTFile(os.path.join(testdir,'compare','neb_int_personalized'))
        self.assertEquals(rname, compare_pr.data)
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

