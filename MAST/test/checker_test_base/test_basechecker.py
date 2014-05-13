"""Tests for Basechecker"""

from MAST.ingredients.checker.basechecker import BaseChecker

import unittest
from unittest import SkipTest
import os
import time
import MAST
import pymatgen
from MAST.utility import dirutil
from MAST.utility import MASTError
from MAST.utility import MASTFile
import shutil
import subprocess
testname="checker_test_base"
print "************************CURDIR:", os.getcwd()
testdir = dirutil.get_test_dir(testname)

class TestBasechecker(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)
        if not os.path.isdir("childdir"):
            os.mkdir("childdir")

    def tearDown(self):
        if os.path.isdir("childdir"):
            shutil.rmtree("childdir")

    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(allowed_keys, **kwargs)

    def test_is_complete(self):
        raise SkipTest
        #self.testclass.is_complete()

    def test_is_started(self):
        raise SkipTest
        #self.testclass.is_started()

    def test_is_ready_to_run(self):
        raise SkipTest
        #self.testclass.is_ready_to_run()

    def test_get_coordinates_only_structure_from_input(self):
        raise SkipTest
        #self.testclass.get_coordinates_only_structure_from_input()

    def test_softlink_a_file(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        mybc.softlink_a_file('childdir','alphatest')
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/alphatest" in listres)
        #self.testclass.softlink_a_file(childpath, filename)
    def test_softlink_a_file_emptyparentfile(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        mybc.softlink_a_file('childdir','emptyfile')
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/emptyfile" in listres)
    def test_softlink_a_file_no_parent_file(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        self.assertRaises(MASTError,mybc.softlink_a_file,'childdir','nofile')
    def test_softlink_a_file_child_file_exists(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        shutil.copy('files/another_alphatest','childdir/alphatest')
        mybc=BaseChecker(allowed,name='files')
        mybc.softlink_a_file('childdir','alphatest')
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertFalse("files/alphatest" in listres)
        self.assertTrue("childdir/alphatest" in myfiles)
        #self.testclass.softlink_a_file(childpath, filename)
    def test_softlink_a_file_child_softlink_exists(self):
        print "CURDIR: ", os.getcwd()
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        os.chdir("childdir")
        linkme=subprocess.Popen("ln -s ../files/another_alphatest alphatest",shell=True)
        linkme.wait()
        os.chdir(testdir)
        mybc=BaseChecker(allowed,name='files')
        mybc.softlink_a_file('childdir','alphatest')
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/another_alphatest" in listres)
        self.assertTrue("childdir/alphatest" in myfiles)
        #self.testclass.softlink_a_file(childpath, filename)
    def test_copy_a_file(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        mybc.copy_a_file('childdir','alphatest','childalpha')
        self.assertTrue(os.path.isfile('childdir/childalpha'))
    def test_copy_a_file_emptyparentfile(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        mybc.copy_a_file('childdir','emptyfile','childalpha')
        self.assertTrue(os.path.isfile('childdir/childalpha'))
    def test_copy_a_file_no_parent_file(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        mybc=BaseChecker(allowed,name='files')
        self.assertRaises(MASTError,mybc.copy_a_file,'childdir','nofile','childalpha')
    def test_copy_a_file_child_file_exists(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        shutil.copy('files/another_alphatest','childdir/alphatest')
        mybc=BaseChecker(allowed,name='files')
        mybc.copy_a_file('childdir','third_alphatest','alphatest')
        compare=MASTFile("files/another_alphatest")
        myfile=MASTFile("childdir/alphatest")
        self.assertEqual(compare.data, myfile.data)
    def test_copy_a_file_child_softlink_exists(self):
        allowed=dict()
        allowed['name']=(str,"","Directory name")
        os.chdir("childdir")
        linkme=subprocess.Popen("ln -s ../files/another_alphatest alphatest",shell=True)
        linkme.wait()
        os.chdir(testdir)
        mybc=BaseChecker(allowed,name='files')
        mybc.copy_a_file('childdir','third_alphatest','alphatest')
        myfiles=dirutil.walkfiles("childdir")
        print myfiles
        #self.assertTrue("childdir/CHGCAR" in myfiles)
        listme=subprocess.Popen("ls -l childdir",stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        listres=listme.communicate()[0]
        listme.wait()
        print listres
        self.assertTrue("files/another_alphatest" in listres)
        self.assertTrue("childdir/alphatest" in myfiles)
        #self.testclass.softlink_a_file(childpath, filename)
