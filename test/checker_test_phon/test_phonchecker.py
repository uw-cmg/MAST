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
        for fname in ["POSCAR","POSCAR_prePHON","INPHON"]:
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
        poscar_pre = MASTFile("files/POSCAR_prePHON")
        poscar_pre.to_file("childdir/POSCAR_prePHON")
        phonposcar = MASTFile("files/POSCAR_for_PHON")
        mypc=PhonChecker(name="childdir")
        mypc._phon_poscar_setup()
        myposcar = MASTFile("childdir/POSCAR")
        self.assertEqual(myposcar.data, phonposcar.data)

    def test__phon_inphon_get_non_mast_keywords(self):
        kdict=dict()
        kdict['mast_program']='phon'
        kdict['LFREE']='True'
        kdict['QA']=11
        mypc=PhonChecker(name="childdir",program_keys = kdict)
        mdict=mypc._phon_inphon_get_non_mast_keywords()
        self.assertEqual(mdict, {'LFREE':'True','QA':11})
        #self.testclass._phon_inphon_get_non_mast_keywords()

    def test__phon_inphon_get_allowed_keywords(self):
        kdict=dict()
        kdict['mast_program']='phon'
        kdict['LFREE']='True'
        kdict['QA']=11
        kdict['IBRION']=2 #This is a VASP, not PHON, keyword
        kdict['nosuchkeyword']="hello"
        mypc=PhonChecker(name="childdir",program_keys = kdict)
        mdict=mypc._phon_inphon_get_non_mast_keywords()
        self.assertEqual(mdict, {'LFREE':'True','QA':11})
        #self.testclass._phon_inphon_get_allowed_keywords(allowedpath)

    def test__phon_inphon_setup(self):
        kdict=dict()
        kdict['mast_program']='phon'
        kdict['TEMPERATURE']='1173'
        kdict['LFREE']='.TRUE.'
        kdict['QA']='11'
        kdict['QB']='11'
        kdict['QC']='11'
        kdict['LSUPER']='.FALSE.'
        kdict['NTYPES']='3'
        kdict['MASS']='55.845 20.000 77.123'
        kdict['IPRINT']='0'
        kdict['ND']='3'
        kdict['IBRION']=2 #This is a VASP, not PHON, keyword
        kdict['nosuchkeyword']="hello"
        mypc=PhonChecker(name="childdir",program_keys = kdict)
        mdict=mypc._phon_inphon_setup()
        compare_inphon = MASTFile("files/INPHON")
        my_inphon = MASTFile("childdir/INPHON")
        self.assertItemsEqual(compare_inphon.data, my_inphon.data)
        poscar_pre = MASTFile("files/POSCAR_prePHON")
        poscar_pre.to_file("childdir/POSCAR_prePHON")
        kdict.pop("NTYPES")
        kdict.pop("MASS")
        mypc=PhonChecker(name="childdir",program_keys = kdict)
        mdict=mypc._phon_inphon_setup()
        compare_inphon = MASTFile("files/INPHON_automass")
        my_inphon = MASTFile("childdir/INPHON")
        self.assertItemsEqual(compare_inphon.data, my_inphon.data)
        poscar_pre = MASTFile("files/POSCAR_prePHON_triple")
        poscar_pre.to_file("childdir/POSCAR_prePHON")
        mypc=PhonChecker(name="childdir",program_keys = kdict)
        mdict=mypc._phon_inphon_setup()
        compare_inphon = MASTFile("files/INPHON_automass_triple")
        my_inphon = MASTFile("childdir/INPHON")
        self.assertItemsEqual(compare_inphon.data, my_inphon.data)

        #self.testclass._phon_inphon_setup()

    def test__phon_inphon_get_masses(self):
        mypc=PhonChecker(name="files")
        [nline,mline]=mypc._phon_inphon_get_masses()
        self.assertEqual(nline,"NTYPES=1")
        self.assertEqual(mline,"MASS=26.9815386 ")
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

