"""Tests for Phonchecker"""

from MAST.ingredients.checker.phonchecker import PhonChecker

from unittest import SkipTest
import os
import time
import MAST
import pymatgen
import unittest
from MAST.utility import dirutil
from MAST.utility import MASTFile

testname="checker_test_phon"
testdir = os.path.join(os.getenv("MAST_INSTALL_PATH"),'test',testname)

class TestPhonChecker(unittest.TestCase):

    def setUp(self):
        os.chdir(testdir)

    def tearDown(self):
        for fname in ["POSCAR","POSCAR_prePHON"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass
    def test___init__(self):
        raise SkipTest
        #self.testclass.__init__(**kwargs)

    def test_is_complete(self):
        mypc=PhonChecker(name="ready")
        self.assertFalse(mypc.is_complete())
        mypc=PhonChecker(name="started")
        self.assertFalse(mypc.is_complete())
        mypc=PhonChecker(name="done_thermo")
        self.assertTrue(mypc.is_complete())
        mypc=PhonChecker(name="done_freq")
        self.assertTrue(mypc.is_complete())

        #self.testclass.is_complete()

    def test_is_ready_to_run(self):
        mypc=PhonChecker(name="notready1")
        self.assertFalse(mypc.is_ready_to_run())
        mypc=PhonChecker(name="notready2")
        self.assertFalse(mypc.is_ready_to_run())
        mypc=PhonChecker(name="notready3")
        self.assertFalse(mypc.is_ready_to_run())
        mypc=PhonChecker(name="notready4")
        self.assertFalse(mypc.is_ready_to_run())
        mypc=PhonChecker(name="ready")
        self.assertTrue(mypc.is_ready_to_run())
        mypc=PhonChecker(name="started") #should these actually return false?
        self.assertTrue(mypc.is_ready_to_run())
        mypc=PhonChecker(name="done_thermo")
        self.assertTrue(mypc.is_ready_to_run())
        mypc=PhonChecker(name="done_freq")
        self.assertTrue(mypc.is_ready_to_run())

    def test__phon_poscar_setup(self):
        phonposcar = MASTFile("files/POSCAR_for_PHON")
        mypc=PhonChecker(name="childdir")
        mypc._phon_poscar_setup()
        myposcar = MASTFile("childdir/POSCAR")
        self.assertEqual(myposcar.data, phonposcar.data)

        
        #self.testclass._phon_poscar_setup()

    def test__phon_inphon_get_non_mast_keywords(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_non_mast_keywords()

    def test__phon_inphon_get_allowed_keywords(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_allowed_keywords(allowedpath)

    def test__phon_inphon_setup(self):
        raise SkipTest
        #self.testclass._phon_inphon_setup()

    def test__phon_inphon_get_masses(self):
        raise SkipTest
        #self.testclass._phon_inphon_get_masses()

    def test__phon_forces_setup(self):
        raise SkipTest
        #self.testclass._phon_forces_setup()

    def test__nosd_my_dynmat(self):
        raise SkipTest
        #self.testclass._nosd_my_dynmat()

    def test__replace_my_displacements(self):
        raise SkipTest
        #self.testclass._replace_my_displacements()

    def test_set_up_program_input(self):
        raise SkipTest
        #self.testclass.set_up_program_input()

    def test_is_started(self):
        raise SkipTest
        #self.testclass.is_started()

