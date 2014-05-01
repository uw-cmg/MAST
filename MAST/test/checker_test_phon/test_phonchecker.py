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
        if not os.path.isdir("childdir"):
            os.mkdir("childdir")

    def tearDown(self):
        for fname in ["POSCAR","POSCAR_prePHON","INPHON","FORCES","DYNMAT","XDATCAR","DYNMAT_mod_1","DYNMAT_mod_2"]:
            try:
                os.remove("childdir/%s" % fname)
            except OSError:
                pass

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
        mystr = pymatgen.io.vaspio.Poscar.from_file("files/POSCAR_prePHON").structure
        myxdat = MASTFile("forces/XDATCAR")
        myxdat.to_file("childdir/XDATCAR")
        mydyn = MASTFile("forces/DYNMAT")
        mydyn.to_file("childdir/DYNMAT")
        mypc=PhonChecker(name="childdir",structure=mystr)
        mypc._phon_forces_setup()
        compare_forces = MASTFile("forces/FORCES_finished")
        myforces = MASTFile("childdir/FORCES")
        self.assertEqual(myforces.data, compare_forces.data)
        #self.testclass._phon_forces_setup()

    def test__nosd_my_dynmat(self):
        mystr = pymatgen.io.vaspio.Poscar.from_file("files/POSCAR_prePHON").structure
        mydyn1 = MASTFile("forces/DYNMAT_mod_1")
        mydyn1.to_file("childdir/DYNMAT_mod_1")
        mypc=PhonChecker(name="childdir",structure=mystr)
        mypc._nosd_my_dynmat()
        mydyn2 = MASTFile("childdir/DYNMAT_mod_2")
        dyn2_compare = MASTFile("forces/DYNMAT_mod_2")
        self.assertEqual(mydyn2.data, dyn2_compare.data)
        #self.testclass._nosd_my_dynmat()

    def test__replace_my_displacements(self):
        mystr = pymatgen.io.vaspio.Poscar.from_file("files/POSCAR_prePHON").structure
        myxdat = MASTFile("forces/XDATCAR")
        myxdat.to_file("childdir/XDATCAR")
        mydyn = MASTFile("forces/DYNMAT")
        mydyn.to_file("childdir/DYNMAT")
        mypc=PhonChecker(name="childdir",structure=mystr)
        mypc._replace_my_displacements()
        mydyn1 = MASTFile("childdir/DYNMAT_mod_1")
        dyn1_compare = MASTFile("forces/DYNMAT_mod_1")
        self.assertEqual(mydyn1.data, dyn1_compare.data)
        #self.testclass._replace_my_displacements()

    def test_set_up_program_input(self):
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
        poscar_pre = MASTFile("files/POSCAR_prePHON")
        poscar_pre.to_file("childdir/POSCAR_prePHON")
        myxdat = MASTFile("forces/XDATCAR")
        myxdat.to_file("childdir/XDATCAR")
        mydyn = MASTFile("forces/DYNMAT")
        mydyn.to_file("childdir/DYNMAT")
        mystr = pymatgen.io.vaspio.Poscar.from_file("files/POSCAR_prePHON").structure
        mypc=PhonChecker(name="childdir",program_keys=kdict,structure=mystr)
        mypc.set_up_program_input()
        forces_compare=MASTFile("forces/FORCES_finished")
        poscar_compare=MASTFile("files/POSCAR_for_PHON")
        inphon_compare=MASTFile("files/INPHON")
        myforces = MASTFile("childdir/FORCES")
        myinphon = MASTFile("childdir/INPHON")
        myposcar = MASTFile("childdir/POSCAR")
        self.assertEqual(myforces.data,forces_compare.data)
        self.assertEqual(myposcar.data,poscar_compare.data)
        self.assertItemsEqual(myinphon.data,inphon_compare.data)

        #self.testclass.set_up_program_input()

    def test_is_started(self):
        mypc=PhonChecker(name="notready1")
        self.assertFalse(mypc.is_started())
        mypc=PhonChecker(name="ready")
        self.assertFalse(mypc.is_started())
        mypc=PhonChecker(name="started")
        self.assertTrue(mypc.is_started())
        mypc=PhonChecker(name="done_freq")
        self.assertTrue(mypc.is_started())
        mypc=PhonChecker(name="done_thermo")
        self.assertTrue(mypc.is_started())

        #self.testclass.is_started()

